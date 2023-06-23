#include "include/SlidingCqt.h"
#include <iostream>
#include <algorithm>

constexpr int OctaveNumber = 9;
constexpr int BinsPerOctave = 12;
constexpr bool Windowing = false;

const int blockSize = 1024;
const double samplerate = 48000.;
const int nBlocks = 100;

std::vector<double> audioInputBlock(blockSize, 0.);
std::vector<double> audioOutputBlock(blockSize, 0.);
std::vector<std::complex<double>> cqtDomainBuffers[OctaveNumber][BinsPerOctave];

Cqt::SlidingCqt<BinsPerOctave, OctaveNumber, Windowing> cqt;

int main(int argc, char const *argv[])
{    
    cqt.init(samplerate, blockSize);

    for (unsigned i_octave = 0u; i_octave < OctaveNumber; i_octave++)
    {
        // The sample rates and block sizes of each downsampled octave can be accessed
        const double octaveRate = cqt.getOctaveSampleRate(i_octave);
        const int octaveSize = cqt.getOctaveBlockSize(i_octave);
        for (unsigned i_tone = 0u; i_tone < BinsPerOctave; i_tone++)
        {
            cqtDomainBuffers[i_octave][i_tone].resize(octaveSize, {0., 0.});
        }
    }

    for (unsigned i_block = 0; i_block < nBlocks; i_block++)
    {
        for (unsigned i = 0; i < blockSize; i++)
        {
            double t = (i+i_block*blockSize)*(1/samplerate);

            // exponential growth from 32.7Hz to 15.804KHz: https://www.desmos.com/calculator/eff93byluq
            double f = std::pow((1+9.4398654e-5), i+i_block*blockSize) + 31.7032;

            // audioInputBlock.at(i) = sin(t*2*M_PI*f);// + (t>1)*sin(t*2*M_PI*32.7032) + (t>1.5)*sin(t*2*M_PI*15804.3);
            audioInputBlock.at(i) = sin(t*2*M_PI*15804.3);
        }

        cqt.inputBlock(audioInputBlock.data(), blockSize);

        for (unsigned i_octave = 0u; i_octave < OctaveNumber; i_octave++)
        {
            const size_t nSamplesOctave = cqt.getSamplesToProcess(i_octave);
            CircularBuffer<std::complex<double>> *octaveCqtBuffer = cqt.getOctaveCqtBuffer(i_octave);
            for (unsigned i_tone = 0u; i_tone < BinsPerOctave; i_tone++)
            {
                octaveCqtBuffer[i_tone].pullBlock(cqtDomainBuffers[i_octave][i_tone].data(), nSamplesOctave);
                for (size_t i_sample = 0u; i_sample < nSamplesOctave; i_sample++)
                {
                    std::complex<double> sample = cqtDomainBuffers[i_octave][i_tone][i_sample];
                    std::cout << sample.real() << "+" << sample.imag() << "j" << " ";
                }
                std::cout << std::endl;
                octaveCqtBuffer[i_tone].pushBlock(cqtDomainBuffers[i_octave][i_tone].data(), nSamplesOctave);
            }
            std::cout << std::endl;
        }
    }
    

    return 0;
}