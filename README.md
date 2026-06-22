# CCSDSViterbi

A high-performance C++ implementation of the CCSDS (Consultative Committee for Space Data Systems) standard Rate 1/2 and Rate 3/4 Convolutional Encoder with Viterbi Decoder, implementing CCSDS-131.0-B-5 specifications.

## Overview

This library provides production-grade implementations of convolutional encoding and maximum-likelihood sequence estimation (Viterbi) decoding for satellite and deep-space communications. The implementation strictly adheres to CCSDS Blue Book standards and includes support for both standard rate 1/2 and punctured rate 3/4 convolutional codes.

## Mathematical Foundation

### Convolutional Codes

Convolutional codes are linear time-invariant codes characterized by three parameters: $(n, k, K)$ where:
- $n$ = number of output bits per encoding cycle
- $k$ = number of input bits per encoding cycle  
- $K$ = constraint length (memory depth + 1)

This implementation uses a **rate $\frac{1}{2}$ convolutional code** with constraint length $K=7$, meaning:
- Input: 1 information bit per cycle
- Output: 2 encoded bits per cycle
- Memory: 6 previous information bits stored in shift register

### Generator Polynomials

The encoder uses two generator polynomials specified in the CCSDS standard:

$$G_1(D) = 171_8 = 1111001_2$$
$$G_2(D) = 133_8 = 1011011_2$$

The output bits are computed via parity-check operations:
$$c_1[n] = \bigoplus_{i} u[n-i] \cdot g_1[i]$$
$$c_2[n] = \bigoplus_{i} u[n-i] \cdot g_2[i]$$

where $\oplus$ denotes XOR and $u[n]$ is the input bit stream.

**Note:** Per CCSDS 131.0-B-5 Section 3.3.1(5), the $G_2$ output is inverted to improve spectral properties. This implementation applies and removes the inversion within the convolutional layer.

### Rate 3/4 Puncturing

Rate 3/4 coding is achieved through **periodic puncturing** of the rate 1/2 code. For every 3 input bits, a puncturing pattern is applied:

$$\text{Pattern} = \begin{bmatrix} 1 & 1 \\ 1 & 0 \\ 0 & 1 \end{bmatrix}$$

This yields 4 output bits for every 3 input bits, achieving rate $\frac{3}{4}$.

### Viterbi Decoding Algorithm

The Viterbi decoder performs **maximum-likelihood sequence detection** on the received (possibly corrupted) encoded bits. The algorithm operates on a trellis structure with:

- **States:** $2^{K-1} = 64$ (representing all possible shift register contents)
- **Transitions:** 2 per state (corresponding to input bits 0 and 1)
- **Branch Metrics:** Hamming distance between received bits and expected code symbols
- **Path Metrics:** Accumulated minimum cost from initial state to current state

The algorithm proceeds in $N$ stages (one per input bit) and produces an estimate of the transmitted bit sequence with near-optimal performance.

## Implementation Details

### Specification

- **Standard:** CCSDS 131.0-B-5 (Recommended Standards for Space Data Systems)
- **Encoder Type:** Convolutional (Rate 1/2 and Rate 3/4 via puncturing)
- **Constraint Length:** $K=7$ (6 memory elements)
- **Decoder Type:** Hard-decision Viterbi
- **Output Inversion:** G2 path inversion per CCSDS spec applied and removed internally

### Key Features

- **Frame Synchronization Compatibility:** Decoded bits are CCSDS-compliant and match the plain ASM (Attached Synchronization Marker) 0x1ACFFC1D for frame sync after decoding
- **Efficient Shift Register:** Compact 7-bit shift register implementation
- **Puncturing Support:** Rotating phase counter for consistent rate 3/4 puncturing
- **G2 Inversion Handling:** Automatic application and removal per CCSDS Blue Book
- **Vector-Based I/O:** Efficient batch processing of bit sequences

### File Structure

- `Viterbi.hpp` - Main header containing `ConvolutionalEncoder` and `ViterbiDecoder` classes
- `Log.h` - Logging utilities (dependency)

## Usage

### Encoding

```cpp
#include "Viterbi.hpp"
using namespace sdr::mdm;

// Create encoder with rate 1/2 configuration
ConvolutionalEncoder encoder;
ConvConfig config;
config.p34 = false;  // false = rate 1/2, true = rate 3/4
encoder.Assemble(config);

// Encode input bits
std::vector<uint8_t> inputBits = {0, 1, 0, 1, 1, 0};
std::vector<uint8_t> encodedBits;
encoder.EncodeBits(&inputBits, &encodedBits);
// For rate 1/2: encodedBits will contain 12 bits
```

### Decoding

```cpp
// Create and initialize decoder
ViterbiDecoder decoder;
ConvConfig config;
config.p34 = false;  // Match encoder configuration
decoder.Assemble(config);

// Decode received bits
std::vector<uint8_t> receivedBits = {...};  // Possibly corrupted encoded bits
std::vector<uint8_t> decodedBits;
decoder.Decode(&receivedBits, &decodedBits);
```

## Performance Characteristics

- **Coding Gain:** ~5 dB at $10^{-5}$ bit error rate (rate 1/2)
- **Memory Requirements:** $O(2^K)$ for trellis structure (512 bytes baseline)
- **Computational Complexity:** $O(N \cdot 2^K)$ where $N$ is message length

## References

- CCSDS 131.0-B-5: *Recommendation for Space Data Systems Standards — Telecommand — Part 1: Deep Space Craft.*
- Viterbi, A., "Error bounds for convolutional codes and an asymptotically optimum decoding algorithm," *IEEE Transactions on Information Theory*, vol. 13, no. 2, pp. 260-269, April 1967.

## Author

J. Enrique Peraza, Trivium Solutions LLC

## License

Refer to the LICENSE file for licensing information.
