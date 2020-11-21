#include "cqt.h"

namespace Math {
	template<typename T>
	T Clamp(T a, T min, T max)
	{
		if (a < min) {
			a = min;
		}
		else if (a > max) {
			a = max;
		}

		return a;
	}
}

namespace CQT {

	std::unordered_map<std::string, std::function<int(float*&, int)>> winMap = {
		{"hann", hanning},
		{"hanning", hanning},
		{"hamm", hamming},
		{"hamming", hamming}
	};

	int hanning(float*& x, int M)
	{
		x = new float[M];
		float a = Math::TWO_PI / ((float)M - 1);
		for (int i = 0; i < M; ++i) {
			x[i] = (0.5f - 0.5f * cos(a * i));
		}
		return M;
	}
	int hamming(float*& x, int M)
	{
		x = new float[M];
		float a = Math::TWO_PI / ((float)M - 1);
		for (int i = 0; i < M; ++i) {
			x[i] = (0.54f - 0.46f * cos(a * i));
		}
		return M;
	}


	std::function<int(float*&, int)> getWindow(std::string name)
	{
		return winMap.at(name);
	}

	/**
	* @fn
	* Constant-Qの中心周波数を計算する
	*
	* @param n_bins	Constant-Qのビン数
	* @param fmin 最小の周波数
	* @param bins_per_octave オクターブあたりのビン数
	* @param tuning
	*
	* @return 配列の大きさ
	*/
	int cqtFrequencies(
		float*& freqs,
		const int n_bins,
		const float fmin,
		const int bins_per_octave,
		const float tuning)
	{
		float correction = pow(2.0f, (tuning / bins_per_octave));
		freqs = new float[n_bins];

		for (int i = 0; i < n_bins; ++i) {
			freqs[i] = correction * fmin * std::pow(2.0f, (i / (float)bins_per_octave));
		}

		return n_bins;
	}

	/**
	* @fn
	* Constant-QTのフィルターを計算する
	*
	* @param freqs	Constant-Qの周波数
	* @param winFunc 窓関数
	* @param nfreq 周波数の個数
	* @param sr サンプリング周波数
	* @param Q Q値
	*
	* @return 配列の大きさ
	*/
	int cqtFilter(
		float**& ksin,
		float**& kcos,
		float*& freqs,
		const std::function<int(float*&, int)> winFunc,
		const int nfreq,
		const int sr,
		const int Q)
	{
		// フィルター用配列を確保
		ksin = new float* [nfreq];
		kcos = new float* [nfreq];

		// 周波数ごとに計算
		for (int k = 0; k < nfreq; ++k) {
			// 周波数
			float freq = freqs[k];
			// サンプル数
			int nsample = int(round((sr * Q) / freq));

			// 配列の確保
			ksin[k] = new float[nsample];
			kcos[k] = new float[nsample];

			// 窓関数を取得
			float* window;
			winFunc(window, nsample);

			// 事前に計算しておく
			float a = Math::TWO_PI * freq / sr;
			// サンプルの数だけ計算する
			for (int i = 0; i < nsample; ++i) {
				float t = a * i;
				ksin[k][i] = std::sin(t) * window[i];
				kcos[k][i] = std::cos(t) * window[i];
			}
			// 窓関数を解放
			delete[] window;
		}
		return nfreq;
	}
	/**
	* @fn
	* オーディオ信号の定数Q変換を計算する。
	*
	* @param y	オーディオデータ
	* @param sr サンプリング周波数
	* @param blocks_per_seconds 一秒あたりの時間解像度
	* @param n_bins 周波数ビン数
	* @param bins_per_octave オクターブあたりのビン数
	* @param fmin 最小の周波数
	* @param Qfactor
	* @param sparsity
	* @param window 窓関数
	* @param tuning
	*
	*/
	std::vector<std::vector<float>> cqt(
		const std::vector<float>& y,
		std::string window,
		const int sr,
		const int n_bins,
		const int bins_per_octave,
		const int samples_per_second,
		float fmin,
		const float Qfactor,
		const float tuning)
	{

		// HOP長を計算
		int hopLength = (int)round((1.0 / samples_per_second) * sr);

		// 最小周波数をチューニングに合わせる
		fmin = fmin * pow(2.0f, (tuning / bins_per_octave));

		// 中心周波数の配列
		float* freqs;
		// 中心周波数を求める
		cqtFrequencies(freqs, n_bins, fmin, bins_per_octave);

		float fratio = 1.0f / bins_per_octave;

		// Q値の計算
		const int Q = int((1.0f / (pow(2.0f, fratio) - 1)) * (Qfactor * fratio));

		// 音声の長さ
		const int ylengths = (int)y.size();

		// フレーム数
		int nframe = ylengths / hopLength;

		// フィルターの配列
		float** ksin;
		float** kcos;
		// フィルターを計算
		cqtFilter(ksin, kcos, freqs, getWindow(window), n_bins, sr, Q);

		// 出力用配列の確保
		std::vector<std::vector<float>> ret(n_bins, std::vector<float>(nframe, 0.0f));

		for (int k = 0; k < n_bins; ++k) {
			// サンプル数
			int nsample = (int)round((sr * Q) / freqs[k]);
			int hsample = nsample / 2;

			for (int f = 0; f < nframe; ++f) {
				// スライス位置を計算
				int istart = f * hopLength - hsample;
				// 音声スライス位置
				int ystart = Math::Clamp(istart, 0, ylengths);
				// 重みスライス位置
				int winstart = Math::Clamp(ystart - istart, 0, nsample);
				int winend = Math::Clamp(ylengths - istart, 0, nsample);
				int winLen = winend - winstart;

				// 合計用
				float sumSin = 0.0f;
				float sumCos = 0.0f;
				// 音声データと重みを掛ける
				for (int i = 0; i < winLen; ++i) {
					sumSin += y[ystart + i] * ksin[k][winstart + i];
					sumCos += y[ystart + i] * kcos[k][winstart + i];
				}

				// サンプル数で割る
				sumSin /= nsample;
				sumCos /= nsample;
				ret[k][f] = std::sqrt(sumCos * sumCos + sumSin * sumSin);
			}
		}


		// 配列の開放
		for (int k = 0; k < n_bins; ++k) {
			delete[] ksin[k];
			delete[] kcos[k];
		}
		delete[] ksin;
		delete[] kcos;
		delete[] freqs;

		return ret;
	}

}