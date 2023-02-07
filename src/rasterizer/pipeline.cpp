#include "pipeline.h"

#include "framebuffer.h"
#include "sample_pattern.h"

#include "../lib/log.h"
#include "../lib/mathlib.h"

#include <iostream>

template< PrimitiveType primitive_type, class Program, uint32_t flags >
void Pipeline< primitive_type, Program, flags >::run(
		std::vector< Vertex > const &vertices,
		typename Program::Parameters const &parameters,
		Framebuffer *framebuffer_) {
	//Framebuffer must be non-null:
	assert(framebuffer_);
	auto &framebuffer = *framebuffer_;

	//A1T7: sample loop
	//TODO: update this function to rasterize to *all* sample locations in the framebuffer.
	// This will probably involve inserting a loop of the form:
	//     std::vector< Vec3 > const &samples = framebuffer.sample_pattern.centers_and_weights;
	//     for (uint32_t s = 0; s < samples.size(); ++s) { ... }
	//  around some subset of the code.
	// You will also need to transform the input and output of the rasterize_* functions to
	//   account for the fact they deal with pixels centered at (0.5,0.5).

	std::vector< ShadedVertex > shaded_vertices;
	shaded_vertices.reserve(vertices.size());

	//--------------------------
	//shade vertices:
	for (auto const &v : vertices) {
		ShadedVertex sv;
		Program::shade_vertex( parameters, v.attributes, &sv.clip_position, &sv.attributes );
		shaded_vertices.emplace_back(sv);
	}

	//--------------------------
	//assemble + clip + homogeneous divide vertices:
	std::vector< ClippedVertex > clipped_vertices;

	//reserve some space to avoid reallocations later:
	if constexpr (primitive_type == PrimitiveType::Lines) {
		//clipping lines can never produce more than one vertex per input vertex:
		clipped_vertices.reserve(shaded_vertices.size());
	} else if constexpr (primitive_type == PrimitiveType::Triangles) {
		//clipping triangles can produce up to 8 vertices per input vertex:
		clipped_vertices.reserve(shaded_vertices.size() * 8);
	}

	//coefficients to map from clip coordinates to framebuffer (i.e., "viewport") coordinates:
	//x: [-1,1] -> [0,width]
	//y: [-1,1] -> [0,height]
	//z: [-1,1] -> [0,1] (OpenGL-style depth range)
	Vec3 const clip_to_fb_scale = Vec3{
		framebuffer.width / 2.0f,
		framebuffer.height / 2.0f,
		0.5f
	};
	Vec3 const clip_to_fb_offset = Vec3{
		0.5f * framebuffer.width,
		0.5f * framebuffer.height,
		0.5f
	};

	//helper used to put output of clipping functions into clipped_vertices:
	auto emit_vertex = [&](ShadedVertex const &sv) {
		ClippedVertex cv;
		float inv_w = 1.0f / sv.clip_position.w;
		cv.fb_position = clip_to_fb_scale * inv_w * sv.clip_position.xyz() + clip_to_fb_offset;
		cv.inv_w = inv_w;
		cv.attributes = sv.attributes;
		clipped_vertices.emplace_back(cv);
	};

	//actually do clipping:
	if constexpr (primitive_type == PrimitiveType::Lines) {
		for (uint32_t i = 0; i + 1 < shaded_vertices.size(); i += 2) {
			clip_line( shaded_vertices[i], shaded_vertices[i+1], emit_vertex );
		}
	} else if constexpr (primitive_type == PrimitiveType::Triangles) {
		for (uint32_t i = 0; i + 2 < shaded_vertices.size(); i += 3) {
			clip_triangle( shaded_vertices[i], shaded_vertices[i+1], shaded_vertices[i+2], emit_vertex );
		}
	} else {
		static_assert( primitive_type == PrimitiveType::Lines, "Unsupported primitive type." );
	}


	//--------------------------
	//rasterize primitives:

	std::vector< Fragment > fragments;

	//helper used to put output of rasterization functions into fragments:
	auto emit_fragment = [&](Fragment const &f) {
		fragments.emplace_back(f);
	};
	//actually do rasterization:
	if constexpr (primitive_type == PrimitiveType::Lines) {
		for (uint32_t i = 0; i + 1 < clipped_vertices.size(); i += 2) {
			rasterize_line( clipped_vertices[i], clipped_vertices[i+1], emit_fragment );
		}
	} else if constexpr (primitive_type == PrimitiveType::Triangles) {
		for (uint32_t i = 0; i + 2 < clipped_vertices.size(); i += 3) {
			rasterize_triangle( clipped_vertices[i], clipped_vertices[i+1], clipped_vertices[i+2], emit_fragment );
		}
	} else {
		static_assert( primitive_type == PrimitiveType::Lines, "Unsupported primitive type." );
	}

	//--------------------------
	//depth test + shade + blend fragments:
	uint32_t out_of_range = 0; //check if rasterization produced fragments outside framebuffer (indicates something is wrong with clipping)
	for (auto const &f : fragments) {

		//fragment location (in pixels):
		int32_t x = (int32_t)std::floor(f.fb_position.x);
		int32_t y = (int32_t)std::floor(f.fb_position.y);

		//if clipping is working properly, this condition shouldn't be needed;
		//however, it prevents crashes while you are working on your clipping functions,
		//so we suggest leaving it in place:
		if (x < 0 || (uint32_t)x >= framebuffer.width || y < 0 || (uint32_t)y >= framebuffer.height) {
			++out_of_range;
			continue;
		}

		//local names that refer to destination sample in framebuffer:
		float &fb_depth = framebuffer.depth_at(x,y,0);
		Spectrum &fb_color = framebuffer.color_at(x,y,0);

		//depth test:
		if constexpr ((flags & PipelineMask_Depth) == Pipeline_Depth_Always) {
			//"Always" means the depth test always passes.
		} else if constexpr ((flags & PipelineMask_Depth) == Pipeline_Depth_Never) {
			//"Never" means the depth test never passes.
			continue; //discard this fragment
		} else if constexpr ((flags & PipelineMask_Depth) == Pipeline_Depth_Less) {
			//"Less" means the depth test passes when the new fragment has depth less than the stored depth.
			//A1T4: Depth_Less
			//TODO: implement depth test! We want to only emit fragments that have a depth less than the stored depth, hence "Depth_Less"
			if (f.fb_position.z > fb_depth){
				continue;
			}
		} else {
			static_assert((flags & PipelineMask_Depth) <= Pipeline_Depth_Always, "Unknown depth test flag.");
		}

		//if depth test passes, and depth writes aren't disabled, write depth to depth buffer:
		if constexpr (!(flags & Pipeline_DepthWriteDisableBit)) {
			fb_depth = f.fb_position.z;
		}

		//shade fragment:
		ShadedFragment sf;
		sf.fb_position = f.fb_position;
		Program::shade_fragment(parameters, f.attributes, f.derivatives, &sf.color, &sf.opacity);

		//write color to framebuffer if color writes aren't disabled:
		if constexpr (!(flags & Pipeline_ColorWriteDisableBit)) {

			//blend fragment:
			if constexpr ((flags & PipelineMask_Blend) == Pipeline_Blend_Replace) {
				fb_color = sf.color;
			} else if constexpr ((flags & PipelineMask_Blend) == Pipeline_Blend_Add) {
				//A1T4: Blend_Add
				//TODO: framebuffer color should have fragment color multiplied by fragment opacity added to it.
				fb_color += sf.color*sf.opacity; 
			} else if constexpr ((flags & PipelineMask_Blend) == Pipeline_Blend_Over) {
				//A1T4: Blend_Over
				//TODO: set framebuffer color to the result of "over" blending (also called "alpha blending") the fragment color over the framebuffer color, using the fragment's opacity
				// 		You may assume that the framebuffer color has its alpha premultiplied already, and you just want to compute the resulting composite color
				fb_color = sf.opacity * sf.color + (1-sf.opacity) * fb_color; //<-- replace this line
			} else {
				static_assert((flags & PipelineMask_Blend) <= Pipeline_Blend_Over, "Unknown blending flag.");
			}
		}
	}

	if (out_of_range > 0) {
		if constexpr (primitive_type == PrimitiveType::Lines) {
			warn("Produced %d fragments outside framebuffer; this indicates something is likely wrong with the clip_line function.", out_of_range);
		} else if constexpr (primitive_type == PrimitiveType::Triangles) {
			warn("Produced %d fragments outside framebuffer; this indicates something is likely wrong with the clip_triangle function.", out_of_range);
		}
	}

	

}

//-------------------------------------------------------------------------
//clipping functions

//helper to interpolate between vertices:
template< PrimitiveType p, class P, uint32_t F >
auto Pipeline< p, P, F >::lerp(ShadedVertex const &a, ShadedVertex const &b, float t) -> ShadedVertex {
	ShadedVertex ret;
	ret.clip_position = (b.clip_position - a.clip_position) * t + a.clip_position;
	for (uint32_t i = 0; i < ret.attributes.size(); ++i) {
		ret.attributes[i] = (b.attributes[i] - a.attributes[i]) * t + a.attributes[i];
	}
	return ret;
}

/*
 * clip_line - clip line to portion with -w <= x,y,z <= w, emit vertices of clipped line (if non-empty)
 *  va, vb: endpoints of line
 *  emit_vertex: call to produce truncated line
 *
 * If clipping shortens the line, attributes of the shortened line should respect the pipeline's interpolation mode.
 * 
 * If no portion of the line remains after clipping, emit_vertex will not be called.
 *
 * The clipped line should have the same direction as the full line.
 *
 */
template< PrimitiveType p, class P, uint32_t flags >
void Pipeline< p, P, flags >::clip_line(
		ShadedVertex const &va, ShadedVertex const &vb,
		std::function< void(ShadedVertex const &) > const &emit_vertex
	) {
	//Determine portion of line over which:
	// pt = (b-a) * t + a
	// -pt.w <= pt.x <= pt.w
	// -pt.w <= pt.y <= pt.w
	// -pt.w <= pt.z <= pt.w

	//... as a range [min_t, max_t]:

	float min_t = 0.0f;
	float max_t = 1.0f;
	
	// want to set range of t for a bunch of equations like:
	//    a.x + t * ba.x <= a.w + t * ba.w
	// so here's a helper:
	auto clip_range = [&min_t, &max_t](float l, float dl, float r, float dr) {
		//restrict range such that:
		//l + t * dl <= r + t * dr
		//re-arranging:
		// l - r <= t * (dr - dl)
		if (dr == dl) {
			//want: l - r <= 0
			if (l - r > 0.0f) {
				//works for none of range, so make range empty:
				min_t = 1.0f; max_t = 0.0f;
			}
		} else if (dr > dl) {
			//since dr - dl is positive:
			//want: (l - r) / (dr - dl) <= t
			min_t = std::max(min_t, (l - r) / (dr - dl));
		} else { //dr < dl
			//since dr - dl is negative:
			//want: (l - r) / (dr - dl) >= t
			max_t = std::min(max_t, (l - r) / (dr - dl));
		}
	};
	
	//local names for clip positions and their difference:
	Vec4 const &a = va.clip_position;
	Vec4 const &b = vb.clip_position;
	Vec4 const ba = b-a;

	// -a.w - t * ba.w <= a.x + t * ba.x <= a.w + t * ba.w
	clip_range(-a.w,-ba.w, a.x, ba.x);
	clip_range( a.x, ba.x, a.w, ba.w);
	// -a.w - t * ba.w <= a.y + t * ba.y <= a.w + t * ba.w
	clip_range(-a.w,-ba.w, a.y, ba.y);
	clip_range( a.y, ba.y, a.w, ba.w);
	// -a.w - t * ba.w <= a.z + t * ba.z <= a.w + t * ba.w
	clip_range(-a.w,-ba.w, a.z, ba.z);
	clip_range( a.z, ba.z, a.w, ba.w);

	if (min_t < max_t) {
		if (min_t == 0.0f) {
			emit_vertex(va);
		} else {
			ShadedVertex out = lerp(va,vb,min_t);
			//don't interpolate attributes if in flat shading mode:
			if constexpr ((flags & PipelineMask_Interp) == Pipeline_Interp_Flat) out.attributes = va.attributes;
			emit_vertex(out);
		}
		if (max_t == 1.0f) {
			emit_vertex(vb);
		} else {
			ShadedVertex out = lerp(va,vb,max_t);
			//don't interpolate attributes if in flat shading mode:
			if constexpr ((flags & PipelineMask_Interp) == Pipeline_Interp_Flat) out.attributes = va.attributes;
			emit_vertex(out);
		}
	}
}


/*
 * clip_triangle - clip triangle to portion with -w <= x,y,z <= w, emit resulting shape as triangles (if non-empty)
 *  va, vb, vc: vertices of triangle
 *  emit_vertex: call to produce clipped triangles (three calls per triangle)
 *
 * If clipping truncates the triangle, attributes of the new vertices should respect the pipeline's interpolation mode.
 * 
 * If no portion of the triangle remains after clipping, emit_vertex will not be called.
 *
 * The clipped triangle(s) should have the same winding order as the full triangle.
 *
 */
template< PrimitiveType p, class P, uint32_t flags >
void Pipeline< p, P, flags >::clip_triangle(
		ShadedVertex const &va, ShadedVertex const &vb, ShadedVertex const &vc,
		std::function< void(ShadedVertex const &) > const &emit_vertex
	) {
	//A1EC: clip_triangle
	//TODO: correct code!
	emit_vertex(va);
	emit_vertex(vb);
	emit_vertex(vc);
}

//-------------------------------------------------------------------------
//rasterization functions

/*
 * rasterize_line:
 * calls emit_fragment( frag ) for every pixel "covered" by the line (va.fb_position.xy, vb.fb_position.xy).
 *
 *    a pixel (x,y) is "covered" by the line if it exits the inscribed diamond:
 * 
 *        (x+0.5,y+1)
 *        /        \
 *    (x,y+0.5)  (x+1,y+0.5)
 *        \        /
 *         (x+0.5,y)
 *
 *    to avoid ambiguity, we consider diamonds to contain their left and bottom points
 *    but not their top and right points. 
 * 
 * 	  since 45 degree lines breaks this rule, our rule in general is to rasterize the line as if its
 *    endpoints va and vb were at va + (e, e^2) and vb + (e, e^2) where no smaller nonzero e produces 
 *    a different rasterization result.
 *
 * for each such diamond, pass Fragment frag to emit_fragment, with:
 *  - frag.fb_position.xy set to the center (x+0.5,y+0.5)
 *  - frag.fb_position.z interpolated linearly between va.fb_position.z and vb.fb_position.z
 *  - frag.attributes set to va.attributes (line will only be used in Interp_Flat mode)
 *  - frag.derivatives set to all (0,0)
 *
 * when interpolating the depth (z) for the fragments, you may use any depth the line takes within the pixel
 * (i.e., you don't need to interpolate to, say, the closest point to the pixel center)
 *
 * If you wish to work in fixed point, check framebuffer.h for useful information about the framebuffer's dimensions.
 *
 */

template< PrimitiveType p, class P, uint32_t flags >
void Pipeline< p, P, flags >::rasterize_line(
		ClippedVertex const &va, ClippedVertex const &vb,
		std::function< void(Fragment const &) > const &emit_fragment
	) {
	if constexpr ((flags & PipelineMask_Interp) != Pipeline_Interp_Flat) {
		assert(0 && "rasterize_line should only be invoked in flat interpolation mode.");
	}
	//A1T2: rasterize_line

	//TODO: Check out the block comment above this function for more information on how to fill in this function!
	// 		The OpenGL specification section 3.5 may also come in handy.



	float x1;
	float y1;
	float x2;
	float y2;
	int type_of_line;
	int ab;
	float dx;
	float dy;
	float slope;
	float b;
	float x_prime;
	float y_prime;
	float z;

	dx = vb.fb_position.x - va.fb_position.x;
	dy = vb.fb_position.y - va.fb_position.y;
	slope = dy/dx;

	auto interp_z = [](float x, float y, ClippedVertex const &va, ClippedVertex const &vb){
		float dist_a;
		float dist_b;
		float sum_dist;

		dist_a = (va.fb_position.y - y)*(va.fb_position.y - y) + (va.fb_position.x - x)*(va.fb_position.x - x);
		dist_b = (vb.fb_position.y - y)*(vb.fb_position.y - y) + (vb.fb_position.x - x)*(vb.fb_position.x - x);
		sum_dist = dist_a+dist_b;

		return va.fb_position.z*(1 - dist_a/sum_dist) + vb.fb_position.z*(1- dist_b/sum_dist);
	};

	if (dx == 0){
		// line is vertical
		type_of_line = 2;
		if (dy > 0){
			ab = 0;
			x1 = va.fb_position.x;
			x2 = vb.fb_position.x;
			y1 = va.fb_position.y;
			y2 = vb.fb_position.y;
		}
		else {
			ab = 1;
			x1 = vb.fb_position.x;
			x2 = va.fb_position.x;
			y1 = vb.fb_position.y;
			y2 = va.fb_position.y;
		}
	} 
	else if (dy == 0){
		// line is horizontal
		type_of_line = 3;
		if (dx > 0){
			ab = 0;
			x1 = va.fb_position.x;
			x2 = vb.fb_position.x;
			y1 = va.fb_position.y;
			y2 = vb.fb_position.y;
		}
		else {
			ab = 1;
			x1 = vb.fb_position.x;
			x2 = va.fb_position.x;
			y1 = vb.fb_position.y;
			y2 = va.fb_position.y;
		}
	}
	else {
		slope = dy/dx;
		if (dy < 0.0f && dx < 0.0f){
			ab = 1;
			x1 = vb.fb_position.x;
			x2 = va.fb_position.x;
			y1 = vb.fb_position.y;
			y2 = va.fb_position.y;
		}
		else {
			ab = 0;
			x1 = va.fb_position.x;
			x2 = vb.fb_position.x;
			y1 = va.fb_position.y;
			y2 = vb.fb_position.y;
		}

		if (std::abs(slope) <= 1){
			// absolute slope greater than 0 and less than 1
			type_of_line = 0;
		}
		else {
			// absolute slope greater than 1
			type_of_line = 1;
		}
	}


	// float dy = vb.fb_position.y - va.fb_position.y;

	// float dx = vb.fb_position.x - va.fb_position.x;

	// float slope = dy/dx;



	// we are always going to draw from x1,y1 to x2,y2

	// calc b for ease of calculations
	b = y1 - slope * x1;
	
	// handle starting edge case, does the point reside top right or bottom right of diamond


	// types of starting edge cases
	float yc = std::floor(y1) + 0.5f;
	float xc = std::floor(x1) + 0.5f;
	// 5. slope of 0
	if (slope == 0.0f) {

		if (ab == 0){
			if (y1 + x1 - yc - xc <= 0.5f && (-y1 + x1 - yc - xc <= -0.5f)){
				Fragment frag;
				frag.attributes = va.attributes;
				frag.derivatives.fill(Vec2(0.0f, 0.0f));
				z = interp_z(xc, yc, va, vb);
				frag.fb_position = Vec3(xc, yc, z);
				emit_fragment(frag);
			}
		}
		else {
			yc = std::floor(y2) + 0.5f;
			xc = std::floor(x2) + 0.5f;
			if (y2 - x2 - yc + xc <= 0.5f && y2 + x2 - xc - yc <= 0.5f){
				Fragment frag;
				frag.attributes = va.attributes;
				frag.derivatives.fill(Vec2(0.0f, 0.0f));
				z = interp_z(xc, yc, va, vb);
				frag.fb_position = Vec3(xc, yc, z);
				emit_fragment(frag);
			}
		}
	}
	// 6. slope of inf
	else if (dx == 0.0f){
		if (ab == 0){
			if (y1 - x1 - yc + xc <= 0.5f && y1 + x1 - yc - xc <= 0.5f){
				Fragment frag;
				frag.attributes = va.attributes;
				frag.derivatives.fill(Vec2(0.0f, 0.0f));
				z = interp_z(xc, yc, va, vb);
				frag.fb_position = Vec3(xc, yc, z);
				emit_fragment(frag);
			}
		}
		else {
			yc = std::floor(y2) + 0.5f;
			xc = std::floor(x2) + 0.5f;
			if (y2 + x2 - yc - xc >= -0.5f && -y2 + x2 - yc - xc <= -0.5f){
				Fragment frag;
				frag.attributes = va.attributes;
				frag.derivatives.fill(Vec2(0.0f, 0.0f));
				z = interp_z(xc, yc, va, vb);
				frag.fb_position = Vec3(xc, yc, z);
				emit_fragment(frag);
			}
		}
	}
	// 1. pos slope going from a to b
	else if (slope > 0 && ab == 0){
		// inside bottom left triangle
		if (y1 + x1 - yc - xc < -0.5f){
			// if endpoint outisde diamond draw this point
			if ((y2 + x2 - yc - xc >= 0.5f) || (y2 - x2 - yc + xc >= 0.5f) || (-y2 + x2 -yc -xc >= -0.5f)){
				Fragment frag;
				frag.attributes = va.attributes;
				frag.derivatives.fill(Vec2(0.0f, 0.0f));
				z = interp_z(xc, yc, va, vb);
				frag.fb_position = Vec3(xc, yc, z);
				emit_fragment(frag);
			}
		}
		// inside bottom right triangle
		else if (-y1 + x1 - yc - xc >= -0.5){
			// can add in check for if we end in diamond
			if (slope >= 1){
				x_prime = (yc - b)/slope;
				if (std::floor(x_prime) == std::floor(x1)){
					Fragment frag;
					frag.attributes = va.attributes;
					frag.derivatives.fill(Vec2(0.0f, 0.0f));
					z = interp_z(xc, yc, va, vb);
					frag.fb_position = Vec3(xc, yc, z);
					emit_fragment(frag);
					
				}
			}
		}
		// inside top left triangle
		else if (y1 - x1 - yc + xc >= 0.5){
			// can add in check for if we end in diamond
			if (slope <= 1){

				y_prime = slope*(xc) + b;
				if (std::floor(y_prime) == std::floor(y1)){
					Fragment frag;
					frag.attributes = va.attributes;
					frag.derivatives.fill(Vec2(0.0f, 0.0f));
					z = interp_z(xc, yc, va, vb);
					frag.fb_position = Vec3(xc, yc, z);
					emit_fragment(frag);
				}
			}
		}
		// we are in the diamond
		else if (y1 + x1 - yc - xc < 0.5){
			// can add in check for if we end within diamond 
			Fragment frag;
			frag.attributes = va.attributes; 
			frag.derivatives.fill(Vec2(0.0f, 0.0f));
			z = interp_z(xc, yc, va, vb);
			frag.fb_position = Vec3(xc, yc, z);
			emit_fragment(frag);
		}
	}
	// 2. pos slope but going from b to a
	else if (slope > 0 && ab == 1){
		yc = std::floor(y2) + 0.5f;
		xc = std::floor(x2) + 0.5f;
		// inside top right triangle
		if (y2 + x2 - yc - xc >= 0.5f){
			// if endpoint outisde diamond draw this point
			if ((y1 + x1 - yc - xc >= 0.5f) || (y1 - x1 - yc + xc >= 0.5f) || (-y1 + x1 -yc -xc) >= -0.5f){
				Fragment frag;
				frag.attributes = va.attributes;
				frag.derivatives.fill(Vec2(0.0f, 0.0f));
				z = interp_z(xc, yc, va, vb);
				frag.fb_position = Vec3(xc, yc, z);
				emit_fragment(frag);
			}
		}
		// inside bottom right triangle
		else if (-y2 + x2 - yc - xc >= -0.5){
			// can add in check for if we end in diamond
			if (slope <= 1){
				x_prime = (yc - b)/slope;
				if (std::floor(x_prime) == std::floor(x2)){
					Fragment frag;
					frag.attributes = va.attributes;
					frag.derivatives.fill(Vec2(0.0f, 0.0f));
					z = interp_z(xc, yc, va, vb);
					frag.fb_position = Vec3(xc, yc, z);
					emit_fragment(frag);
					
				}
			}
		}
		// inside top left triangle
		else if (y2 - x2 - yc + xc >= 0.5){
			// can add in check for if we end in diamond
			if (slope >= 1){

				y_prime = slope*(xc) + b;
				if (std::floor(y_prime) == std::floor(y2)){
					Fragment frag;
					frag.attributes = va.attributes;
					frag.derivatives.fill(Vec2(0.0f, 0.0f));
					z = interp_z(xc, yc, va, vb);
					frag.fb_position = Vec3(xc, yc, z);
					emit_fragment(frag);
				}
			}
		}
		// we are in the diamond
		else if (y2 + x2 - yc - xc >= -0.5){
			// can add in check for if we end within diamond 
			Fragment frag;
			frag.attributes = va.attributes;
			frag.derivatives.fill(Vec2(0.0f, 0.0f));
			z = interp_z(xc, yc, va, vb);
			frag.fb_position = Vec3(xc, yc, z);
			emit_fragment(frag);
		}
	}
	// 3. negative slope going from a to b
	else if (slope < 0 && ab == 0){
		// inside top left triangle
		if (y1 - x1 - yc + xc >= 0.5f){
			// if endpoint outisde diamond draw this point (could add check for inside the diamond)
			Fragment frag;
			frag.attributes = va.attributes;
			frag.derivatives.fill(Vec2(0.0f, 0.0f));
			z = interp_z(xc, yc, va, vb);
			frag.fb_position = Vec3(xc, yc, z);
			emit_fragment(frag);
			
		}
		// inside bottom left triangle
		else if (y1 + x1 - yc - xc <= -0.5){
			// can add in check for if we end in diamond
			if (std::abs(slope) <= 1){
				x_prime = (yc - b)/slope;
				if (std::floor(x_prime) == std::floor(x1)){
					Fragment frag;
					frag.attributes = va.attributes;
					frag.derivatives.fill(Vec2(0.0f, 0.0f));
					z = interp_z(xc, yc, va, vb);
					frag.fb_position = Vec3(xc, yc, z);
					emit_fragment(frag);
					
				}
			}
		}
		// inside top right triangle
		else if (y1 - x1 - yc + xc >= 0.5){
			// can add in check for if we end in diamond
			if (std::abs(slope) >= 1){

				y_prime = slope*(xc) + b;
				if (std::floor(y_prime) == std::floor(y1)){
					Fragment frag;
					frag.attributes = va.attributes;
					frag.derivatives.fill(Vec2(0.0f, 0.0f));
					z = interp_z(xc, yc, va, vb);
					frag.fb_position = Vec3(xc, yc, z);
					emit_fragment(frag);
				}
			}
		}
		// we are in the diamond
		else if (-y1 + x1 - yc - xc <= -0.5){
			// can add in check for if we end within diamond 
			Fragment frag;
			frag.attributes = va.attributes;
			frag.derivatives.fill(Vec2(0.0f, 0.0f));
			z = interp_z(xc, yc, va, vb);
			frag.fb_position = Vec3(xc, yc, z);
			emit_fragment(frag);
		}
	}
	// 4. negative slope going from b to a
	else if (slope < 0 && ab == 1){
		yc = std::floor(y2) + 0.5f;
		xc = std::floor(x2) + 0.5f;
		// inside bottom right triangle
		if (-y2 + x2 - yc - xc >= -0.5f){
			// if endpoint outisde diamond draw this point (could add check for inside the diamond)
			Fragment frag;
			frag.attributes = va.attributes;
			frag.derivatives.fill(Vec2(0.0f, 0.0f));
			z = interp_z(xc, yc, va, vb);
			frag.fb_position = Vec3(xc, yc, z);
			emit_fragment(frag);
			
		}
		// inside bottom left triangle
		else if (y2 + x2 - yc - xc <= -0.5){
			// can add in check for if we end in diamond
			if (std::abs(slope) >= 1){
				x_prime = (yc - b)/slope;
				if (std::floor(x_prime) == std::floor(x2)){
					Fragment frag;
					frag.attributes = va.attributes;
					frag.derivatives.fill(Vec2(0.0f, 0.0f));
					z = interp_z(xc, yc, va, vb);
					frag.fb_position = Vec3(xc, yc, z);
					emit_fragment(frag);
					
				}
			}
		}
		// inside top right triangle
		else if (y2 - x2 - yc + xc >= 0.5){
			// can add in check for if we end in diamond
			if (std::abs(slope) <= 1){

				y_prime = slope*(xc) + b;
				if (std::floor(y_prime) == std::floor(y2)){
					Fragment frag;
					frag.attributes = va.attributes;
					frag.derivatives.fill(Vec2(0.0f, 0.0f));
					z = interp_z(xc, yc, va, vb);
					frag.fb_position = Vec3(xc, yc, z);
					emit_fragment(frag);
				}
			}
		}
		// we are in the diamond
		else if (y2 - x2 - yc + xc <= 0.5){
			// can add in check for if we end within diamond 
			Fragment frag;
			frag.attributes = va.attributes;
			frag.derivatives.fill(Vec2(0.0f, 0.0f));
			z = interp_z(xc, yc, va, vb);
			frag.fb_position = Vec3(xc, yc, z);
			emit_fragment(frag);
		}
	}





	// for loop to emit intermediate fragments (multiple cases)

	// horizontal case
	if (type_of_line == 3) {
		y_prime = std::floor(y1);
		for (x_prime = std::floor(x1) + 1.0f; x_prime <= std::floor(x2) - 1.0f; x_prime++){
			Fragment frag;
			frag.attributes = va.attributes;
			frag.derivatives.fill(Vec2(0.0f, 0.0f));
			z = interp_z(x_prime + 0.5f, y_prime + 0.5f, va, vb);
			frag.fb_position = Vec3(x_prime + 0.5f, y_prime + 0.5f, z);
			emit_fragment(frag);
		}

	}
	// vertical case
	else if (type_of_line == 2){
		x_prime = std::floor(x1);

		for (y_prime = std::floor(y1) + 1.0f; y_prime <= std::floor(y2) - 1.0f; y_prime++){
			Fragment frag;
			frag.attributes = va.attributes;
			frag.derivatives.fill(Vec2(0.0f, 0.0f));
			z = interp_z(x_prime + 0.5f, y_prime + 0.5f, va, vb);
			frag.fb_position = Vec3(x_prime + 0.5f, y_prime + 0.5f, z);
			emit_fragment(frag);
		}
	}
	// pos abs slope greater than 1
	else if (type_of_line == 1){
		for (y_prime = std::floor(y1) + 1.0f; y_prime <= std::floor(y2) - 1.0f; y_prime++){
			x_prime = (y_prime + 0.5f - b)/slope;
			Fragment frag;  
			frag.attributes = va.attributes;
			frag.derivatives.fill(Vec2(0.0f, 0.0f));
			z = interp_z(std::floor(x_prime) + 0.5f, std::floor(y_prime) + 0.5f, va, vb);
			frag.fb_position = Vec3(std::floor(x_prime) + 0.5f, y_prime + 0.5f, z);
			emit_fragment(frag);
		}
	}
	// pos abs slope less than 1 case
	else if (type_of_line == 0){
		// may not need to floor x2
		for (x_prime = std::floor(x1) + 1.0f; x_prime <= std::floor(x2) - 1.0f; x_prime++){
			y_prime = slope*(x_prime + 0.5f) + b;
			Fragment frag;
			frag.attributes = va.attributes;
			frag.derivatives.fill(Vec2(0.0f, 0.0f));
			z = interp_z(std::floor(x_prime) + 0.5f, std::floor(y_prime) + 0.5f, va, vb);
			frag.fb_position = Vec3(x_prime + 0.5f, std::floor(y_prime) + 0.5f, z);
			emit_fragment(frag);
		}
	}

	// handle end edge case, does the point reside top right or bottom right of diamond
	// 1. slope is 0 
	yc = std::floor(y2) + 0.5f;
	xc = std::floor(x2) + 0.5f;
	if (slope == 0){
		if (ab == 0){
			if ((y2 + x2 - yc - xc >= 0.5f) || (-y2 + x2 - yc - xc <= -0.5f)){
				Fragment frag;
				frag.attributes = va.attributes;
				frag.derivatives.fill(Vec2(0.0f, 0.0f));
				z = interp_z(xc, yc, va, vb);
				frag.fb_position = Vec3(xc, yc, z);
				emit_fragment(frag);
			}
		}
		else {
			yc = std::floor(y1) + 0.5f;
			xc = std::floor(x1) + 0.5f;
			if (y1 - x1 - yc + xc >= 0.5 || y1 + x1 - yc - xc <= -0.5){
				Fragment frag;
				frag.attributes = va.attributes;
				frag.derivatives.fill(Vec2(0.0f, 0.0f));
				z = interp_z(xc, yc, va, vb);
				frag.fb_position = Vec3(xc, yc, z);
				emit_fragment(frag);
			}
		}
	}
	// slope is inf
	else if (dx == 0){
		if (ab == 0){
			if (y2 - x2 - yc + xc >= 0.5f || y2 + x2 - yc - xc >= 0.5f){
				Fragment frag;
				frag.attributes = va.attributes;
				frag.derivatives.fill(Vec2(0.0f, 0.0f));
				z = interp_z(xc, yc, va, vb);
				frag.fb_position = Vec3(xc, yc, z);
				emit_fragment(frag);
			}
		}
		else {
			yc = std::floor(y1) + 0.5f;
			xc = std::floor(x1) + 0.5f;
			if (y1 + x1 - yc - xc <= -0.5f || -y1 + x1 - yc - xc >= -0.5f){
				Fragment frag;
				frag.attributes = va.attributes;
				frag.derivatives.fill(Vec2(0.0f, 0.0f));
				z = interp_z(xc, yc, va, vb);
				frag.fb_position = Vec3(xc, yc, z);
				emit_fragment(frag);
			}
		}
	}
	// slope is pos and a to b
	else if (slope > 0 && ab == 0){
		// inside top right triangle
		if (y2 + x2 - yc - xc >= 0.5f){
			// if startpoint outisde diamond draw this point ( could be checked )
			Fragment frag;
			frag.attributes = va.attributes;
			frag.derivatives.fill(Vec2(0.0f, 0.0f));
			z = interp_z(xc, yc, va, vb);
			frag.fb_position = Vec3(xc, yc, z);
			emit_fragment(frag);
		}
		// inside bottom right triangle
		else if (-y2 + x2 - yc - xc >= -0.5){
			// can add in check for if we end in diamond
			if (slope <= 1){
				x_prime = (yc - b)/slope;
				if (std::floor(x_prime) == std::floor(x2)){
					Fragment frag;
					frag.attributes = va.attributes;
					frag.derivatives.fill(Vec2(0.0f, 0.0f));
					z = interp_z(xc, yc, va, vb);
					frag.fb_position = Vec3(xc, yc, z);
					emit_fragment(frag);
					
				}
			}
		}
		// inside top left triangle
		else if (y2 - x2 - yc + xc >= 0.5){
			// can add in check for if we end in diamond
			if (slope >= 1){
				y_prime = slope*(xc) + b;
				if (std::floor(y_prime) == std::floor(y2)){
					Fragment frag;
					frag.attributes = va.attributes;
					frag.derivatives.fill(Vec2(0.0f, 0.0f));
					z = interp_z(xc, yc, va, vb);
					frag.fb_position = Vec3(xc, yc, z);
					emit_fragment(frag);
				}
			}
		}
	}
	// slope is pos and b to a
	else if (slope > 0 && ab == 1){
		yc = std::floor(y1) + 0.5f;
		xc = std::floor(x1) + 0.5f;
		// inside bottom left triangle
		if (y1 + x1 - yc - xc <= -0.5f){
			// if startpoint outisde diamond draw this point ( could be checked )
			Fragment frag;
			frag.attributes = va.attributes;
			frag.derivatives.fill(Vec2(0.0f, 0.0f));
			z = interp_z(xc, yc, va, vb);
			frag.fb_position = Vec3(xc, yc, z);
			emit_fragment(frag);
		}
		// inside bottom right triangle
		else if (-y1+ x1 - yc - xc >= -0.5){
			// can add in check for if we end in diamond
			if (slope >= 1){
				x_prime = (yc - b)/slope;
				if (std::floor(x_prime) == std::floor(x1)){
					Fragment frag;
					frag.attributes = va.attributes;
					frag.derivatives.fill(Vec2(0.0f, 0.0f));
					z = interp_z(xc, yc, va, vb);
					frag.fb_position = Vec3(xc, yc, z);
					emit_fragment(frag);
					
				}
			}
		}
		// inside top left triangle
		else if (y1 - x1 - yc + xc >= 0.5){
			// can add in check for if we end in diamond
			if (slope <= 1){

				y_prime = slope*(xc) + b;
				if (std::floor(y_prime) == std::floor(y2)){
					Fragment frag;
					frag.attributes = va.attributes;
					frag.derivatives.fill(Vec2(0.0f, 0.0f));
					z = interp_z(xc, yc, va, vb);
					frag.fb_position = Vec3(xc, yc, z);
					emit_fragment(frag);
				}
			}
		}
	}
	// slope negative a to b
	else if (slope < 0 && ab == 0){
		// inside bottom right triangle
		if (-y1 + x1 - yc - xc >= -0.5f){
			// if endpoint outisde diamond draw this point (could add check for inside the diamond)
			Fragment frag;
			frag.attributes = va.attributes;
			frag.derivatives.fill(Vec2(0.0f, 0.0f));
			z = interp_z(xc, yc, va, vb);
			frag.fb_position = Vec3(xc, yc, z);
			emit_fragment(frag);
			
		}
		// inside bottom left triangle
		else if (y2 + x2 - yc - xc <= -0.5){
			// can add in check for if we end in diamond
			if (std::abs(slope) >= 1){
				x_prime = (yc - b)/slope;
				if (std::floor(x_prime) == std::floor(x2)){
					Fragment frag;
					frag.attributes = va.attributes;
					frag.derivatives.fill(Vec2(0.0f, 0.0f));
					z = interp_z(xc, yc, va, vb);
					frag.fb_position = Vec3(xc, yc, z);
					emit_fragment(frag);
					
				}
			}
		}
		// inside top right triangle
		else if (y2 - x2 - yc + xc >= 0.5){
			// can add in check for if we end in diamond
			if (std::abs(slope) <= 1){

				y_prime = slope*(xc) + b;
				if (std::floor(y_prime) == std::floor(y2)){
					Fragment frag;
					frag.attributes = va.attributes;
					frag.derivatives.fill(Vec2(0.0f, 0.0f));
					z = interp_z(xc, yc, va, vb);
					frag.fb_position = Vec3(xc, yc, z);
					emit_fragment(frag);
				}
			}
		}
	}
	// slope negative b to a
	else if (slope < 0 && ab == 1){
		yc = std::floor(y1) + 0.5f;
		xc = std::floor(x1) + 0.5f;
		// inside top left triangle
		if (y1 - x1 - yc + xc >= 0.5f){
			// if endpoint startpoint diamond draw this point (could add check for inside the diamond)
			Fragment frag;
			frag.attributes = va.attributes;
			frag.derivatives.fill(Vec2(0.0f, 0.0f));
			z = interp_z(xc, yc, va, vb);
			frag.fb_position = Vec3(xc, yc, z);
			emit_fragment(frag);
			
		}
		// inside bottom left triangle
		else if (y1 + x1 - yc - xc <= -0.5){
			// can add in check for if we end in diamond
			if (std::abs(slope) <= 1){
				x_prime = (yc - b)/slope;
				if (std::floor(x_prime) == std::floor(x1)){
					Fragment frag;
					frag.attributes = va.attributes;
					frag.derivatives.fill(Vec2(0.0f, 0.0f));
					z = interp_z(xc, yc, va, vb);
					frag.fb_position = Vec3(xc, yc, z);
					emit_fragment(frag);
					
				}
			}
		}
		// inside top right triangle
		else if (y1 - x1 - yc + xc >= 0.5){
			// can add in check for if we end in diamond
			if (std::abs(slope) >= 1){
				y_prime = slope*(xc) + b;
				if (std::floor(y_prime) == std::floor(y1)){
					Fragment frag;
					frag.attributes = va.attributes;
					frag.derivatives.fill(Vec2(0.0f, 0.0f));
					z = interp_z(xc, yc, va, vb);
					frag.fb_position = Vec3(xc, yc, z);
					emit_fragment(frag);
				}
			}
		}
	}


}


/*
 *
 * rasterize_triangle(a,b,c,emit) calls 'emit(frag)' at every location
 *  (x+0.5,y+0.5) (where x,y are integers) covered by triangle (a,b,c).
 *
 * The emitted fragment should have:
 * - frag.fb_position.xy = (x+0.5, y+0.5)
 * - frag.fb_position.z = linearly interpolated fb_position.z from a,b,c (NOTE: does not depend on Interp mode!)
 * - frag.attributes = depends on Interp_* flag in flags:
 *   - if Interp_Flat: copy from va.attributes
 *   - if Interp_Screen: interpolate as if (a,b,c) is a 2D triangle flat on the screen
 *   - if Interp_Correct: use perspective-correct interpolation
 * - frag.derivatives = derivatives w.r.t. fb_position.x and fb_position.y of the first frag.derivatives.size() attributes.
 *
 * Notes on derivatives:
 *  The derivatives are partial derivatives w.r.t. screen locations. That is:
 *    derivatives[i].x = d/d(fb_position.x) attributes[i]
 *    derivatives[i].y = d/d(fb_position.y) attributes[i]
 *  You may compute these derivatives analytically or numerically.
 *
 *  See section 8.12.1 "Derivative Functions" of the GLSL 4.20 specification for some inspiration. (*HOWEVER*, the spec is solving a harder problem, and also nothing in the spec is binding on your implementation)
 *
 *  One approach is to rasterize blocks of four fragments and use forward and backward differences to compute derivatives.
 *  To assist you in this approach, keep in mind that the framebuffer size is *guaranteed* to be even. (see framebuffer.h)
 *
 * Notes on coverage:
 *  If two triangles are on opposite sides of the same edge, and a
 *  fragment center lies on that edge, rasterize_triangle should
 *  make sure that exactly one of the triangles emits that fragment.
 *  (Otherwise, speckles or cracks can appear in the final render.)
 * 
 *  For degenerate (co-linear) triangles, you may consider them to not be on any side of an edge.
 * 	Thus, even if two degnerate triangles share an edge that contains a fragment center, you don't need to emit it.
 *  You will not lose points for doing something reasonable when handling this case
 *
 *  This is pretty tricky to get exactly right!
 *
 */
template< PrimitiveType p, class P, uint32_t flags >
void Pipeline< p, P, flags >::rasterize_triangle(
		ClippedVertex const &va, ClippedVertex const &vb, ClippedVertex const &vc,
		std::function< void(Fragment const &) > const &emit_fragment
	) {
	//NOTE: it is okay to restructure this function to allow these tasks to use the
	// same code paths. Be aware, however, that all of them need to remain working!
	// (e.g., if you break Flat while implementing Correct, you won't get points
	//  for Flat.)

	float min_x = std::floor(std::min(va.fb_position.x, std::min(vb.fb_position.x, vc.fb_position.x)));
	float max_x = std::ceil(std::max(va.fb_position.x, std::max(vb.fb_position.x, vc.fb_position.x)));
	float min_y = std::floor(std::min(va.fb_position.y, std::min(vb.fb_position.y, vc.fb_position.y)));
	float max_y = std::ceil(std::max(va.fb_position.y, std::max(vb.fb_position.y, vc.fb_position.y)));
	

	float det1;
	float det2;
	float det3;
	float bary_a;
	float bary_b;
	float bary_c;
	bool signs;
	int edge;
	float z; 
	float det_sum;
	
	auto determine_sign = [](float cx, float cy, ClippedVertex const &p1, ClippedVertex const &p2) {
		return (cx - p2.fb_position.x) * (p1.fb_position.y - p2.fb_position.y) - (p1.fb_position.x - p2.fb_position.x) * (cy - p2.fb_position.y);
	};

	auto matching_signs3 = [](float f1, float f2, float f3){
		return (f1 > 0.0f && f2 > 0.0f && f3 > 0.0f) || (f1 < 0.0f && f2 < 0.0f && f3 < 0.0f);
	};

	auto edge_case = [](bool match, float d1, float d2, float d3) {
		// det1 associated with c
		// det2 associated with a
		// det3 associated with b
		if (!match){
			if(!d1 && !d2){
				// falls on b
				return 4;
			}
			else if (!d1 && !d3){
				// falls on a
				return 5;
			}
			else if (!d2 && !d3){
				// falls on c
				return 6;
			}
			else if (!d1 && ((d2 > 0.0f && d3 > 0.0f) || (d2 < 0.0f && d3 < 0.0f))){
				return 1;
			}
			else if (!d2 && ((d1 > 0.0f && d3 > 0.0f) || (d1 < 0.0f && d3 < 0.0f))){
				return 2;
			}
			else if (!d3 && ((d2 > 0.0f && d1 > 0.0f) || (d2 < 0.0f && d1 < 0.0f))){
				return 3;
			}
		}
		return 0;
	};

	auto point_on_edge = [](ClippedVertex const &non_edge, ClippedVertex const &edge1, ClippedVertex const &edge2){
		float dy = edge2.fb_position.y - edge1.fb_position.y;
		float dx = edge2.fb_position.x - edge1.fb_position.x;

		if (dx == 0){
			if (non_edge.fb_position.x > edge1.fb_position.x){
				return true;
			}
		}
		else if (dy == 0){
			if (non_edge.fb_position.y < edge1.fb_position.y){
				return true;
			}
		}
		
		float slope = dy/dx;
		float b = edge1.fb_position.y - slope*(edge1.fb_position.x);
		if (slope > 0){
			if (non_edge.fb_position.y < slope*non_edge.fb_position.x + b) {
				return true;
			}
		}
		else if (slope < 0) {
			if (non_edge.fb_position.y > slope*non_edge.fb_position.x + b){
				return true;
			}
		}

		return false;
	};

	auto interpolate_attr = [&](float ptx, float pty, ClippedVertex const &va, ClippedVertex const &vb, ClippedVertex const &vc){
		std::array< float, FA > ret_attr = va.attributes;
		float d1 = determine_sign(ptx, pty, va, vb);
		float d2 = determine_sign(ptx, pty, vb, vc);
		float d3 = determine_sign(ptx, pty, vc, va);
		float d_sum = d1 + d2 + d3;
		for(int x = 0; x < 5; x++){
			ret_attr[x] = (d1/d_sum) * vc.attributes[x] + (d2/d_sum) * va.attributes[x] + (d3/d_sum) * vb.attributes[x];
		}
		return ret_attr;
	};

	auto interpolate_attr2 = [&](float ptx, float pty, std::array< float, FA > const &va_attr, 
	std::array< float, FA > const &vb_attr, std::array< float, FA > const &vc_attr){
		std::array< float, FA > ret_attr = va.attributes;
		float d1 = determine_sign(ptx, pty, va, vb);
		float d2 = determine_sign(ptx, pty, vb, vc);
		float d3 = determine_sign(ptx, pty, vc, va);
		float d_sum = d1 + d2 + d3;
		for(int x = 0; x < 5; x++){
			ret_attr[x] = (d1/d_sum) * vc_attr[x] + (d2/d_sum) * va_attr[x] + (d3/d_sum) * vb_attr[x];
		}
		return ret_attr;
	};

	auto interpolate_z = [&](float ptx, float pty, ClippedVertex const &va, ClippedVertex const &vb, ClippedVertex const &vc){
		float d1 = determine_sign(ptx, pty, va, vb);
		float d2 = determine_sign(ptx, pty, vb, vc);
		float d3 = determine_sign(ptx, pty, vc, va);
		float d_sum = d1 + d2 + d3;
		return (d1/d_sum) * vc.inv_w + (d2/d_sum) * va.inv_w + (d3/d_sum) * vb.inv_w;
	};

	if constexpr ((flags & PipelineMask_Interp) == Pipeline_Interp_Flat) {
		//A1T3: flat triangles
		//TODO: rasterize triangle (see block comment above this function).

		// iterate from bottom left to top right emitting all fragments that are inside the triangle
		for (float curr_x = min_x + 0.5f; curr_x <= max_x + 0.5; curr_x++){
			for (float curr_y = min_y + 0.5f; curr_y <= max_y + 0.5; curr_y++){
				
				det1 = determine_sign(curr_x, curr_y, va, vb);
				det2 = determine_sign(curr_x, curr_y, vb, vc);
				det3 = determine_sign(curr_x, curr_y, vc, va);
				
				signs = matching_signs3(det1, det2, det3);

				edge = edge_case(signs, det1, det2, det3);

				det_sum = det1 + det2 + det3;
				z = (det1/det_sum)*(vc.fb_position.z) + (det2/det_sum)*(va.fb_position.z) + (det3/det_sum)*(vb.fb_position.z);

				if (edge != 0){
					bool emit = false;
					// todo fill in for each type of edge case
					if (edge == 1){
						// c is the point not on the edge that our test point falls on
						emit = point_on_edge(vc, va, vb);
					}
					else if (edge == 2){
						// a is the point not on the edge that our test point falls on
						emit = point_on_edge(va, vc, vb);
					}
					else if (edge == 3){
						// b is the point not on the edge that our test point falls on 
						emit = point_on_edge(vb, va, vc);
					}
					else if (edge == 4){
						// falls on vertex b
						emit = point_on_edge(vc, va, vb) && point_on_edge(va, vc, vb);
					}
					else if (edge == 5){
						// falls on vertex a
						emit = point_on_edge(vb, va, vc) && point_on_edge(vc, va, vb);
					}
					else if (edge == 6){
						// falls on vertex c
						emit = point_on_edge(va, vc, vb) && point_on_edge(vb, vc, va);
					}

					if (emit){
						Fragment frag;
						frag.attributes = va.attributes;
						frag.derivatives.fill(Vec2(0.0f, 0.0f));
						frag.fb_position = Vec3(curr_x, curr_y, z);
						emit_fragment(frag);
					}
				}
				else if (signs == 1){
					Fragment frag;
					frag.attributes = va.attributes;
					frag.derivatives.fill(Vec2(0.0f,0.0f));
					frag.fb_position = Vec3(curr_x, curr_y, z);
					emit_fragment(frag);
				}
				
			}
		}

	} else if constexpr ((flags & PipelineMask_Interp) == Pipeline_Interp_Screen) {
		//A1T5: screen-space smooth triangles
		//TODO: rasterize triangle (see block comment above this function).

		// iterate from bottom left to top right emitting all fragments that are inside the triangle
		for (float curr_x = min_x + 0.5f; curr_x <= max_x + 0.5; curr_x++){
			for (float curr_y = min_y + 0.5f; curr_y <= max_y + 0.5; curr_y++){
				
				det1 = determine_sign(curr_x, curr_y, va, vb);
				det2 = determine_sign(curr_x, curr_y, vb, vc);
				det3 = determine_sign(curr_x, curr_y, vc, va);
				
				signs = matching_signs3(det1, det2, det3);

				edge = edge_case(signs, det1, det2, det3);

				det_sum = det1 + det2 + det3;

				bary_a = det2/det_sum;
				bary_b = det3/det_sum;
				bary_c = det1/det_sum;
					

				z = (bary_c)*(vc.fb_position.z) + (bary_a)*(va.fb_position.z) + (bary_b)*(vb.fb_position.z);

				if (edge != 0){
					bool emit = false;
					// todo fill in for each type of edge case
					if (edge == 1){
						// c is the point not on the edge that our test point falls on
						emit = point_on_edge(vc, va, vb);
					}
					else if (edge == 2){
						// a is the point not on the edge that our test point falls on
						emit = point_on_edge(va, vc, vb);
					}
					else if (edge == 3){
						// b is the point not on the edge that our test point falls on 
						emit = point_on_edge(vb, va, vc);
					}
					else if (edge == 4){
						// falls on vertex b
						emit = point_on_edge(vc, va, vb) && point_on_edge(va, vc, vb);
					}
					else if (edge == 5){
						// falls on vertex a
						emit = point_on_edge(vb, va, vc) && point_on_edge(vc, va, vb);
					}
					else if (edge == 6){
						// falls on vertex c
						emit = point_on_edge(va, vc, vb) && point_on_edge(vb, vc, va);
					}

					if (emit){
						Fragment frag;
						frag.fb_position = Vec3(curr_x, curr_y, z);
						frag.attributes = interpolate_attr(curr_x, curr_y, va, vb, vc);
						std::array< float, FA > temp_x = interpolate_attr(curr_x + 1, curr_y, va, vb, vc);
						std::array< float, FA > temp_y = interpolate_attr(curr_x, curr_y + 1, va, vb, vc);
						frag.derivatives.fill(Vec2(0.0f, 0.0f));
						Vec2 deriv_1 = Vec2(temp_x[0] - frag.attributes[0], temp_x[1] - frag.attributes[1]);
						Vec2 deriv_2 = Vec2(temp_y[0] - frag.attributes[0], temp_y[1] - frag.attributes[1]);
						frag.derivatives[0] = deriv_1;
						frag.derivatives[1] = deriv_2;
						emit_fragment(frag);
					}
				}
				else if (signs == 1){
					Fragment frag;
					frag.fb_position = Vec3(curr_x, curr_y, z);
					frag.attributes = interpolate_attr(curr_x, curr_y, va, vb, vc);
					std::array< float, FA > temp_x = interpolate_attr(curr_x + 1, curr_y, va, vb, vc);
					std::array< float, FA > temp_y = interpolate_attr(curr_x, curr_y + 1, va, vb, vc);
					frag.derivatives.fill(Vec2(0.0f, 0.0f)); 
					Vec2 deriv_1 = Vec2(temp_x[0] - frag.attributes[0], temp_y[0] - frag.attributes[0]);
					Vec2 deriv_2 = Vec2(temp_x[1] - frag.attributes[1], temp_y[1] - frag.attributes[1]);
					frag.derivatives[0] = deriv_1;
					frag.derivatives[1] = deriv_2;
					emit_fragment(frag);
				}
				
			}
		}

		



	} else if constexpr ((flags & PipelineMask_Interp) == Pipeline_Interp_Correct) {
		//A1T5: perspective correct triangles
		//TODO: rasterize triangle (block comment above this function).

		

		for (float curr_x = min_x + 0.5f; curr_x <= max_x + 0.5; curr_x++){
			for (float curr_y = min_y + 0.5f; curr_y <= max_y + 0.5; curr_y++){
				
				det1 = determine_sign(curr_x, curr_y, va, vb);
				det2 = determine_sign(curr_x, curr_y, vb, vc);
				det3 = determine_sign(curr_x, curr_y, vc, va);
				
				signs = matching_signs3(det1, det2, det3);

				edge = edge_case(signs, det1, det2, det3);

				det_sum = det1 + det2 + det3;
				z = (det1/det_sum)*(vc.fb_position.z) + (det2/det_sum)*(va.fb_position.z) + (det3/det_sum)*(vb.fb_position.z);

				if (edge != 0){
					bool emit = false;
					// todo fill in for each type of edge case
					if (edge == 1){
						// c is the point not on the edge that our test point falls on
						emit = point_on_edge(vc, va, vb);
					}
					else if (edge == 2){
						// a is the point not on the edge that our test point falls on
						emit = point_on_edge(va, vc, vb);
					}
					else if (edge == 3){
						// b is the point not on the edge that our test point falls on 
						emit = point_on_edge(vb, va, vc);
					}
					else if (edge == 4){
						// falls on vertex b
						emit = point_on_edge(vc, va, vb) && point_on_edge(va, vc, vb);
					}
					else if (edge == 5){
						// falls on vertex a
						emit = point_on_edge(vb, va, vc) && point_on_edge(vc, va, vb);
					}
					else if (edge == 6){
						// falls on vertex c
						emit = point_on_edge(va, vc, vb) && point_on_edge(vb, vc, va);
					}

					if (emit){
						Fragment frag;

						// multiply inv_w by attr
						std::array< float, FA > attr_a = va.attributes;
						std::array< float, FA > attr_b = va.attributes;
						std::array< float, FA > attr_c = va.attributes;
						std::array< float, FA > frag_attr;
						float frag_z;

						for (int i = 0; i < 5; i++){
							attr_a[i] = va.inv_w * va.attributes[i];
							attr_b[i] = vb.inv_w * vb.attributes[i];
							attr_c[i] = vc.inv_w * vc.attributes[i];
						}

						frag_attr = interpolate_attr2(curr_x, curr_y, attr_a, attr_b, attr_c);
						frag_z = interpolate_z(curr_x, curr_y, va, vb, vc);
						
						for(int i = 0; i < 5; i++){
							frag_attr[i] /= frag_z;
						}

						frag.attributes = frag_attr;

						std::array< float, FA > temp_x = interpolate_attr2(curr_x + 1, curr_y, attr_a, attr_b, attr_c);
						float tmp_x_z = interpolate_z(curr_x + 1, curr_y, va, vb, vc);
						std::array< float, FA > temp_y = interpolate_attr2(curr_x, curr_y + 1, attr_a, attr_b, attr_c);
						float tmp_y_z = interpolate_z(curr_x, curr_y + 1, va, vb, vc);

						for(int i = 0; i < 5; i++){
							temp_x[i] /= tmp_x_z;
							temp_y[i] /= tmp_y_z;
						}

						frag.derivatives.fill(Vec2(0.0f, 0.0f));
						Vec2 deriv1 = Vec2(temp_x[0] - frag.attributes[0], temp_y[0] - frag.attributes[0]);
						Vec2 deriv2 = Vec2(temp_x[1] - frag.attributes[1], temp_y[1] - frag.attributes[1]);
						frag.derivatives[0] = deriv1;
						frag.derivatives[1] = deriv2;
						frag.fb_position = Vec3(curr_x, curr_y, z);
						emit_fragment(frag);
					}
				}
				else if (signs == 1){
					Fragment frag;

						// multiply inv_w by attr
						std::array< float, FA > attr_a = va.attributes;
						std::array< float, FA > attr_b = va.attributes;
						std::array< float, FA > attr_c = va.attributes;
						std::array< float, FA > frag_attr;
						float frag_z;

						for (int i = 0; i < 5; i++){
							attr_a[i] = va.inv_w * va.attributes[i];
							attr_b[i] = vb.inv_w * vb.attributes[i];
							attr_c[i] = vc.inv_w * vc.attributes[i];
						}

						frag_attr = interpolate_attr2(curr_x, curr_y, attr_a, attr_b, attr_c);
						frag_z = interpolate_z(curr_x, curr_y, va, vb, vc);
						
						for(int i = 0; i < 5; i++){
							frag_attr[i] /= frag_z;
						}

						frag.attributes = frag_attr;

						std::array< float, FA > temp_x = interpolate_attr2(curr_x + 1, curr_y, attr_a, attr_b, attr_c);
						float tmp_x_z = interpolate_z(curr_x + 1, curr_y, va, vb, vc);
						std::array< float, FA > temp_y = interpolate_attr2(curr_x, curr_y + 1, attr_a, attr_b, attr_c);
						float tmp_y_z = interpolate_z(curr_x, curr_y + 1, va, vb, vc);

						for(int i = 0; i < 5; i++){
							temp_x[i] /= tmp_x_z;
							temp_y[i] /= tmp_y_z;
						}

						frag.derivatives.fill(Vec2(0.0f, 0.0f));
						Vec2 deriv1 = Vec2(temp_x[0] - frag.attributes[0], temp_y[0] - frag.attributes[0]);
						Vec2 deriv2 = Vec2(temp_x[1] - frag.attributes[1], temp_y[1] - frag.attributes[1]);
						frag.derivatives[0] = deriv1;
						frag.derivatives[1] = deriv2;
						frag.fb_position = Vec3(curr_x, curr_y, z);
						emit_fragment(frag);
				}
				
			}
		}

	}
}


//-------------------------------------------------------------------------
//compile instantiations for all programs and blending and testing types:

#include "programs.h"

template struct Pipeline< PrimitiveType::Lines, Programs::Lambertian, Pipeline_Blend_Replace | Pipeline_Depth_Less | Pipeline_Interp_Flat >;
template struct Pipeline< PrimitiveType::Lines, Programs::Lambertian, Pipeline_Blend_Replace | Pipeline_Depth_Always | Pipeline_Interp_Flat >;
template struct Pipeline< PrimitiveType::Triangles, Programs::Lambertian, Pipeline_Blend_Replace | Pipeline_Depth_Less | Pipeline_Interp_Flat >;
template struct Pipeline< PrimitiveType::Triangles, Programs::Lambertian, Pipeline_Blend_Replace | Pipeline_Depth_Less | Pipeline_Interp_Screen >;
template struct Pipeline< PrimitiveType::Triangles, Programs::Lambertian, Pipeline_Blend_Replace | Pipeline_Depth_Less | Pipeline_Interp_Correct >;
