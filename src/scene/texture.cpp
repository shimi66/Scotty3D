
#include "texture.h"

#include <iostream>

namespace Textures {


Spectrum sample_nearest(HDR_Image const &image, Vec2 uv) {
	//clamp texture coordinates, convert to [0,w]x[0,h] pixel space:
	float x = image.w * std::clamp(uv.x, 0.0f, 1.0f);
	float y = image.h * std::clamp(uv.y, 0.0f, 1.0f);

	//the pixel with the nearest center is the pixel that contains (x,y):
	int32_t ix = int32_t(std::floor(x));
	int32_t iy = int32_t(std::floor(y));

	//texture coordinates of (1,1) map to (w,h), and need to be reduced:
	ix = std::min(ix, int32_t(image.w) - 1);
	iy = std::min(iy, int32_t(image.h) - 1);

	return image.at(ix, iy);
}

Spectrum sample_bilinear(HDR_Image const &image, Vec2 uv) {
	//A1T6: sample_bilinear
	//TODO: implement bilinear sampling strategy on texture 'image'
	// std::cout << "uv " << uv << std::endl;
	// std::cout << "image w and h " << image.w << " " << image.h  << std::endl;
	
	float x = image.w * std::clamp(uv.x, 0.0f, 1.0f);
	float y = image.h * std::clamp(uv.y, 0.0f, 1.0f);

	// std::cout << "x, y " << x << " " << y << std::endl;

	int32_t ix = int32_t(std::floor(x - 0.5));
	int32_t iy = int32_t(std::floor(y - 0.5));

	ix = std::min(ix, int32_t(image.w) - 1);
	iy = std::min(iy, int32_t(image.h) - 1);
	// std::cout << "ix and iy " << ix << " " << iy << std::endl;

	float s = x - ix - 0.5f;
	float t = y - iy - 0.5f;
	// std::cout << "s " << s << std::endl;
	// std::cout << "t " << t << std::endl;

	// std::cout << "possible value??? " << image.at(ix+1, iy+1) << std::endl;

	// use current box if i am on top or left edge (4 cases)
	if (ix + 1 == (int) image.w && iy + 1 == (int) image.h){
		return image.at(ix, iy);
	}
	else if (ix + 1 == (int) image.w) {
		return (1-t)*image.at(ix, iy) + t*image.at(ix, iy+1);
	}
	else if (iy + 1 == (int) image.h){
		return (1-s)*image.at(ix, iy) + s*image.at(ix+1, iy);
	}
	else{
		return (1-t)*((1-s)*image.at(ix, iy) + s*image.at(ix+1, iy)) + t*((1-s)*image.at(ix, iy+1) + s*image.at(ix+1, iy+1));
	}

}


Spectrum sample_trilinear(HDR_Image const &base, std::vector< HDR_Image > const &levels, Vec2 uv, float lod) {
	//A1T6: sample_trilinear
	//TODO: implement trilinear sampling strategy on using mip-map 'levels'

	// std::cout << "base w and h " << base.w << " " << base.h  << std::endl;
	// std::cout << "uv " << uv << std::endl;
	// std::cout << "lod " << lod << std::endl;

	// std::cout << "x, y " << x << " " << y << std::endl;

	int32_t iz = int32_t(std::floor(lod));
	// std::cout << "ix and iy iz " << ix << " " << iy << " " << iz << std::endl;

	// call bilinear each level to get h0, h1 (see notes)
	// interpolate between h0, h1
	Spectrum h0;

	if (lod == std::floor(lod)){
		if (lod == 0){ 
			return sample_bilinear(base, uv);
		}
		else {
			return sample_bilinear(levels[(int) lod-1], uv);
		}
	}
	else {
		if (iz == 0){ 
			h0 = sample_bilinear(base, uv);
		}
		else {
			h0 = sample_bilinear(levels[iz-1], uv);
		}
		Spectrum h1 = sample_bilinear(levels[iz], uv);

		// std::cout << "h0 and h1 " << h0 << " " << h1 << std::endl;

		float coeff = lod - std::floor(lod);

		Spectrum output = (1 - coeff)*h0 + coeff*h1;

		// 0.75 =  lod  iz = 0

		return output;
	}
	 
}

/*
 * generate_mipmap- generate mipmap levels from a base image.
 *  base: the base image
 *  levels: pointer to vector of levels to fill (must not be null)
 *
 * generates a stack of levels [1,n] of sizes w_i, h_i, where:
 *   w_i = max(1, floor(w_{i-1})/2)
 *   h_i = max(1, floor(h_{i-1})/2)
 *  with:
 *   w_0 = base.w
 *   h_0 = base.h
 *  and n is the smalles n such that w_n = h_n = 1
 *
 * each level should be calculated by downsampling a blurred version
 * of the previous level to remove high-frequency detail.
 *
 */
void generate_mipmap(HDR_Image const &base, std::vector< HDR_Image > *levels_) {
	assert(levels_);
	auto &levels = *levels_;


	{ // allocate sublevels sufficient to scale base image all the way to 1x1:
		int32_t num_levels = static_cast<int32_t>(std::log2(std::max(base.w, base.h)));
		assert(num_levels >= 0);

		levels.clear();
		levels.reserve(num_levels);

		uint32_t width = base.w;
		uint32_t height = base.h;
		for (int32_t i = 0; i < num_levels; ++i) {
			assert(!(width == 1 && height == 1)); //would have stopped before this if num_levels was computed correctly

			width = std::max(1u, width / 2u);
			height = std::max(1u, height / 2u);

			levels.emplace_back(width, height);
		}
		assert(width == 1 && height == 1);
		assert(levels.size() == uint32_t(num_levels));
	}

	//now fill in the levels using a helper:
	//downsample:
	// fill in dst to represent the low-frequency component of src
	auto downsample = [](HDR_Image const &src, HDR_Image &dst) {
		//dst is half the size of src in each dimension:
		assert(std::max(1u, src.w / 2u) == dst.w);
		assert(std::max(1u, src.h / 2u) == dst.h);

		//A1T6: generate
		std::cout << std::endl;
		//TODO: Write code to fill the levels of the mipmap hierarchy by downsampling
		for (int width = 0; width < (int) dst.w; width++){
			for (int height = 0; height < (int) dst.h; height++){
				Spectrum tmp = src.at(2*width, 2*height);
				tmp += src.at(2*width+1, 2*height);
				tmp += src.at(2*width, 2*height+1);
				tmp += src.at(2*width+1, 2*height+1);
				tmp.r /= 4;
				tmp.g /= 4;
				tmp.b /= 4;
				dst.at(width, height) = tmp;
				std::cout << "dst.at " << dst.at(width, height) << std::endl;
			}
			std::cout << std::endl;
		}

		//Be aware that the alignment of the samples in dst and src will be different depending on whether the image is even or odd.

	};

	std::cout << "Regenerating mipmap (" << levels.size() << " levels): [" << base.w << "x" << base.h << "]";
	std::cout.flush();
	for (uint32_t i = 0; i < levels.size(); ++i) {
		HDR_Image const &src = (i == 0 ? base : levels[i-1]);
		HDR_Image &dst = levels[i];
		std::cout << " -> [" << dst.w << "x" << dst.h << "]"; std::cout.flush();

		downsample(src, dst);
	}
	std::cout << std::endl;
	
}

Image::Image(Sampler sampler_, HDR_Image const &image_) {
	sampler = sampler_;
	image = image_.copy();
	update_mipmap();
}

Spectrum Image::evaluate(Vec2 uv, float lod) const {
	if (sampler == Sampler::nearest) {
		// std::cout << "in nearest " << std::endl;
		return sample_nearest(image, uv);
	} else if (sampler == Sampler::bilinear) {
		// std::cout << "in bilinear " << std::endl;
		return sample_bilinear(image, uv);
	} else {
		// std::cout << "in trilinear " << std::endl;
		return sample_trilinear(image, levels, uv, lod);
	}
}

void Image::update_mipmap() {
	if (sampler == Sampler::trilinear) {
		generate_mipmap(image, &levels);
	} else {
		levels.clear();
	}
}

GL::Tex2D Image::to_gl() const {
	return image.to_gl(1.0f);
}

void Image::make_valid() {
	update_mipmap();
}

Spectrum Constant::evaluate(Vec2 uv, float lod) const {
	return color * scale;
}

} // namespace Textures
bool operator!=(const Textures::Constant& a, const Textures::Constant& b) {
	return a.color != b.color || a.scale != b.scale;
}

bool operator!=(const Textures::Image& a, const Textures::Image& b) {
	return a.image != b.image;
}

bool operator!=(const Texture& a, const Texture& b) {
	if (a.texture.index() != b.texture.index()) return false;
	return std::visit(
		[&](const auto& data) { return data != std::get<std::decay_t<decltype(data)>>(b.texture); },
		a.texture);
}
