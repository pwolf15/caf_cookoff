#include <iostream>
#include <string>
#include <complex>
#include <fstream>
#include <vector>
#include <cstring>

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

std::vector<std::complex<float>> apply_fdoa(std::vector<std::complex<float>> ray, float fdoa, float sample_rate)
{
  std::complex<float> precache(0,2*3.14159265*fdoa/sample_rate);
  std::vector<std::complex<float>> new_ray(ray.size());
  for (size_t i = 0; i < new_ray.size(); ++i)
  {
    new_ray[i] = ray[i] * std::exp(precache*std::complex<float>(i,i));
  }

  return new_ray;
}

std::vector<float> amb_surf(
  std::vector<std::complex<float>>& needle,
  std::vector<std::complex<float>>& haystack, 
  std::vector<float>& freqs_hz,
  float sample_rate)
{
  uint32_t len_needle = needle.size();
  uint32_t len_haystack = haystack.size();
  uint32_t len_freqs = freqs_hz.size();
  std::vector<float> surf(len_freqs*len_needle);
  for (size_t i = 0; i < len_freqs; ++i)
  {
    std::vector<std::complex<float>> shifted = apply_fdoa(needle, freqs_hz[i], sample_rate);
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