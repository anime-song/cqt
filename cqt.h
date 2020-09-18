#pragma once
#include <cmath>
#include <functional>
#include <vector>
#include <unordered_map>

namespace Math {
	constexpr float PI = 3.14159265358979323846f;
	constexpr float TWO_PI = PI * 2.0f;

	template<typename T> T Clamp(T a, T min, T max);

}

namespace CQT {

	int hanning(float*& x, int M);
	int hamming(float*& x, int M);

	static std::unordered_map<const char*, std::function<int(float*&, int)>> winMap;
	std::function<int(float*&, int)> getWindow(const char* name);

	int cqtFrequencies(
		float*& freqs,
		const int n_bins,
		const float fmin,
		const int bins_per_octave = 12,
		const float tuning = 0.0f);

	int cqtFilter(
		float**& ksin,
		float**& kcos,
		float*& freqs,
		const std::function<int(float*&, int)> winFunc,
		const int nfreq,
		const int sr,
		const int Q);

	std::vector<std::vector<float>> cqt(
		const std::vector<float>& y,
		const char* window = "hanning",
		const int sr = 22050,
		const int n_bins = 12 * 3 * 7,
		const int bins_per_octave = 12 * 3,
		const int bins_per_second = 40,
		float fmin = 32.703f,
		const float Qfactor = 20.0f,
		const float tuning = 0.0f);
}
