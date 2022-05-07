#include <iostream>
#include <string>
#include <complex>
#include <complex.h>
#include <fstream>
#include <vector>
#include <cstring>
#include <math.h>
#include "fftw3.h"

void xcorr(std::vector<std::complex<float>>& signala, 
  std::vector<std::complex<float>> signalb, 
  std::vector<std::complex<float>> result, int N)
{
    std::vector<std::complex<float>> signala_ext(2*N-1), signalb_ext(2*N-1), out_shifted(2*N-1);
    std::vector<std::complex<float>> outa(2*N-1), outb(2*N-1), out(2*N-1);

    fftw_plan pa = fftw_plan_dft_1d(2 * N - 1, reinterpret_cast<fftw_complex*>(signala_ext.data()), reinterpret_cast<fftw_complex*>(outa.data()), FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan pb = fftw_plan_dft_1d(2 * N - 1, reinterpret_cast<fftw_complex*>(signalb_ext.data()), reinterpret_cast<fftw_complex*>(outb.data()), FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan px = fftw_plan_dft_1d(2 * N - 1, reinterpret_cast<fftw_complex*>(out.data()), reinterpret_cast<fftw_complex*>(result.data()), FFTW_BACKWARD, FFTW_ESTIMATE);

    //zeropadding
    for (int i = N - 1; i <= N + (N - 1); ++i)
    {
      signala_ext[i] = signala[i - (N - 1)];
    }
    for (int i = N; i <= N + (N - 1); ++i)
    {
      signalb_ext[i] = signalb[i - N];
    }

    fftw_execute(pa);
    fftw_execute(pb);

    std::complex<float> scale = 1.0/(2 * N -1);
    for (int i = 0; i < 2 * N - 1; i++)
        out[i] = outa[i] * conj(outb[i]) * scale;

    fftw_execute(px);

    fftw_destroy_plan(pa);
    fftw_destroy_plan(pb);
    fftw_destroy_plan(px);

    fftw_cleanup();

    return;
}

std::vector<std::complex<float>> read_complex(std::string filename)
{
  std::ifstream file;
  file.open(filename, std::ios::binary);

  std::vector<unsigned char> buffer(std::istreambuf_iterator<char>(file), {});
  std::vector<std::complex<float>> data(buffer.size() / sizeof(std::complex<float>));
  std::cout << buffer.size() << std::endl;

  std::memcpy((char*)data.data(), buffer.data(), buffer.size());

  return data;
};

std::vector<std::complex<float>> apply_fdoa(std::vector<std::complex<float>>& ray, float fdoa, float sample_rate)
{
  std::complex<float> precache(0,2*3.141592653589*fdoa/sample_rate);
  std::vector<std::complex<float>> new_ray(ray.size());
  for (size_t i = 0; i < 10; ++i)
  {
    new_ray[i] = ray[i] * std::exp(precache*std::complex<float>(i,i));
  }

  return new_ray;
}

std::vector<float> xcor(std::vector<std::complex<float>>& apple, std::vector<std::complex<float>>& banana)
{
  std::vector<std::complex<float>> corr(2 * apple.size() - 1);

  xcorr(apple, banana, corr, apple.size());

  std::vector<float> amplitudes;
  for (int i = 0; i < corr.size(); ++i)
  {
    amplitudes.push_back(sqrt(pow(corr[i].real(), 2) + pow(corr[i].imag(), 2)));
  }

  return amplitudes;
}

std::vector<std::vector<float>> amb_surf(
  std::vector<std::complex<float>>& needle,
  std::vector<std::complex<float>>& haystack, 
  std::vector<float>& freqs_hz,
  float sample_rate)
{
  uint32_t len_needle = needle.size();
  uint32_t len_haystack = haystack.size();
  uint32_t len_freqs = freqs_hz.size();
  std::vector<std::vector<float>> surf(len_freqs);
  for (size_t i = 0; i < 10; ++i)
  {
    std::vector<std::complex<float>> shifted = apply_fdoa(needle, freqs_hz[i], sample_rate);
    std::cout << shifted[i] << std::endl;
    surf[i] = xcor(shifted, haystack);
  }
  

  return surf;
}

int main()
{

  std::string data_dir = "../../data";
  std::string needle_filename = "chirp_4_raw.c64";
  std::string haystack_filename = "chirp_4_T+70samp_F+82.89Hz.c64";
  std::cout << haystack_filename << std::endl;

  std::cout << "CAF result" << std::endl;

  std::vector<std::complex<float>> data1 = read_complex(data_dir + "/" + needle_filename);
  std::cout << "Num samples data 1: " << data1.size() << std::endl;

  std::vector<std::complex<float>> data2 = read_complex(data_dir + "/" + haystack_filename);
  data2.resize(data1.size());
  std::cout << "Num samples data 2: " << data2.size() << std::endl;

  float sample_rate = 48e3;
  std::vector<float> freq_offsets;
  float c = -100;
  while (c < 100)
  {
    freq_offsets.push_back(c);
    c += 0.5;
  }

  std::cout << "freq offsets: " << freq_offsets.size() << std::endl;
  amb_surf(data1, data2, freq_offsets, sample_rate);

  return 0;
}