 /* 
* *
* * Filename: ConvolutionalEncoder.hpp
* *
* * Description:
* *   CCSDS-standard Rate 1/2 Convolutional Encoder, K=7, Polynomials 171(octal),133(octal).
* *  - Constraint length K=7 (6 memory elements). Rate 3/4 via puncturing over 3 input bits P=[[1.1.1],[1.0.1]]
* * 
* *
* * Author:
* *   JEP, J. Enrique Peraza
* *
* *
*/
///NOTE: G2 OUTPUT INVERSION (CCSDS 131.0-B-5 sec 3.3.1(5)) is applied and removed
// within the convolutional layer. Post-decode information bits are CCSDS-plain
// so any consumer that frame-syncs on DECODED BITS matches the plain ASM 0x1ACFFC1D and must not
// re-apply he inversion. A consumer that ever correlates ENCODED symbols (frame sync ahead of Viterbi)
// must use the G2-inverted encoded ASM as its reference. This layer guarantees the former;
// it does not produce the latter. 
#pragma once
#include <algorithm>
#include <array>
#include <cstdio>
#include <cstdint>
#include <vector>
#include <array>
#include "Log.h"
namespace sdr::mdm
{
  struct ConvConfig
  {
    bool p34{false};               // True for rate 3/4 via puncturing; false for rate 1/2
    ConvConfig (void)=default;     // Default constructor
  };
  class ConvolutionalEncoder
  {
    public:
      ConvolutionalEncoder (void)
      {
        Reset();                    // Initialize shift register
      }
      ~ConvolutionalEncoder (void)=default;
      inline void Reset (void) { z=0; ph=0; } // Reset shift register and config
      inline void Assemble (const ConvConfig& ccfg)
      {                               // ~~~~~~~~~~ Assemble ~~~~~~~~~~ //
        this->ccfg=ccfg;              // Store configuration
      }                               // ~~~~~~~~~~ Assemble ~~~~~~~~~~ //
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      // Parity-check bit generation for one input byte (8 bits)
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      inline uint8_t Parity (uint32_t bs) // Bit sequence to parity check
      {                               // ~~~~~~~~~~ Parity ~~~~~~~~~~ //                                
        bs^=bs>>16u;                  // Fold upper bits
        bs^=bs>>8u;                   // Fold again
        bs^=bs>>4u;                   // Fold again
        bs^=bs>>2u;                   // Fold again
        bs^=bs>>1u;                   // Fold again
        return static_cast<uint8_t>(bs&0x1u);// Return parity bit
      }                                 // ~~~~~~~~~~ Parity ~~~~~~~~~~ //
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      // Encode bit sequence
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      inline void Encode (
        uint8_t u,                      // Input bit to encode (0 or 1)
        std::vector<uint8_t>* const o)  // Output encoded bits
      {                                 // ~~~~~~~~~~ Encode ~~~~~~~~~~ //
        //Deb(" [BEGIN] Convolutional Encode: Input Bit=%d",u);
        z=((z<<1)|(u&0x01))&0x7F;       // Update shift register with input bit
        //Deb(" Convolutional Encode: Shift Reg=%02X",z);
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        // Generators (octal): 171=0b1111001, 133=0b1011011 (reversed for shift register)
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        uint32_t g1=0b1111001u;         // Generator 1
        uint32_t g2=0b1011011u;         // Generator 2
        uint8_t p1=Parity(z&g1);        // Compute parity bit 1
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        // CCSDS 131-0-B-5 sec 3.3.1(5): symbol inveersion on the G2 output path.
        // Inverted HERE at the source of so rate 1/2 and all rate 3/4 puncture
        // phases emit the inverted C2 consistently
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        uint8_t p2=Parity(z&g2);        // Compute parity bit 2 (G2=133 octal), raw
        if (!ccfg.p34)                  // No puncturing?
          p2=static_cast<uint8_t>(p2^1u);// Invert G2 for the basic rate 1/2 only
        //Deb(" Convolutional Encode: Input Bit=%d Shift Reg=%02X Parity Bits=[%d,%d]",u,z,p1,p2);
        if (!ccfg.p34)                  // Rate 1/2?
        {                               // Yes, output both parity bits
          o->push_back(p1);             // Output parity bit 1
          o->push_back(p2);             // Output parity bit 2
        }                               // Done rate 1/2
        else                            // Else its is punctured rate 3/4
        {                               //
          // ~~~~~~~~~~~~~~~~~~~~~~~~~~ //
          // Rate 3/4 puncture over groups of 3 inputs -> pattern indexes
          // Implemented using a rotating phase counter mod 3
          // ~~~~~~~~~~~~~~~~~~~~~~~~~~ //
          if (ph==0)                    // Phase 0: output both parity bits
          {                             //
            o->push_back(p1);           // Output parity bit 1
            o->push_back(p2);           // Output parity bit 2
          }                             //
          else if (ph==1)               // Phase 1: output only parity bit 1
          {                             //
            o->push_back(p2);           // Output parity bit 2 C2
          }                             //
          else                          // Phase 2: output only parity bit 2
          {                             //
            o->push_back(p1);           // Output parity bit 1
          }                             //
          ph=(ph+1)%3;                  // Advance puncturing phase
          //Deb(" [END] Convolutional Encode: Puncturing Phase=%d",ph);
        }                               // Done rate 3/4
      }                                 // ~~~~~~~~~~ Encode ~~~~~~~~~~ //
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      // Encode a block of input bits
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      inline void EncodeBits (
        const std::vector<uint8_t>* in, // Input bits
        std::vector<uint8_t>* const o)  // Output encoded bits
      {                                 // ~~~~~~~~~~ EncodeBits ~~~~~~~~~~ //
        if (in==nullptr||o==nullptr)    // Bad args?
          return;                       // Nothing to do
        o->clear();                     // Clear output buffer
        o->reserve(in->size()*2);       // Reserve space (max rate 1/2)
        for (const auto& b:(*in))       // For each input bit
          Encode(b,o);                  // Encode bit
      }                                 // ~~~~~~~~~~ EncodeBits ~~~~~~~~~~ //
    private:
      uint32_t z{0};                    // Shift register state
      ConvConfig ccfg;                  // Configuration parameters
      int ph{0};                        // Puncturing phase counter for encoder (0,1,2)
  };
  class ViterbiDecoder
  {
    public:
      ViterbiDecoder (void)
      {                                 // ~~~~~~~~~~ Constructor ~~~~~~~~~~ //
        Assemble(ConvConfig{});         // Default config
      }                                 // ~~~~~~~~~~ Constructor ~~~~~~~~~~ //
      ~ViterbiDecoder (void)
      {
        //Log("ViterbiDecoder destructor called");
      }
      inline void Assemble (const ConvConfig& conf)
      {                                 // ~~~~~~~~~~ Assemble ~~~~~~~~~~ //
        ccfg=conf;                      // Store configuration
        trellis.BuildTrellis();         // Build trellis structure
      }                                 // ~~~~~~~~~~ Assemble ~~~~~~~~~~ //
      inline void Reset (void)
      {                                 // ~~~~~~~~~~ Reset ~~~~~~~~~~ //
        for (int s=0;s<64;++s)
          pm[s]=INF;                    // Initialize path metrics to infinity
        pm[0]=0;                        // Start at state 0
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        // Puncturing phase aligns with encoder (0-> p1,p2; 1 -> p1; 2 -> p2)
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        ph=0;                           // Reset puncturing phase
        steps=0;                        // Reset step counter
        prev.clear();                   // Clear survivor paths
        bit.clear();                    // Clear input bits
        dmarg.clear();                  // Clear SOVA margins
        prevc.clear();                  // Clear SOVA competitor predecessors
        bitc.clear();                   // Clear SOVA competitor input bits
      }                                 // ~~~~~~~~~~ Reset ~~~~~~~~~~ //
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      // Decode a stream of hard bits (0/1). Output recovered input bits.
      // For rate 1/2: ch.size() must be even -> one input bit per 2 channel bits (what the encoder produced)
      // For rate 3/2: uses puncturing phases; every 4 channel bits correspond to 3 input bits.
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      inline void DecodeBits (
        const std::vector<uint8_t>* ch, // Input channel bits (hard 0/1)
        std::vector<uint8_t>* const o)  // Output decoded input bits
      {                                 // ~~~~~~~~~~ DecodeBits ~~~~~~~~~~ //
        if (ch==nullptr||o==nullptr||ch->empty())
          return;
        o->clear();                     // Clear output buffer
        ////Log("[BEGIN] Viterbi DecodeBits: input channel bits=%zu",ch->size());
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        // Prepare survivor path storage
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        prev.reserve(prev.size()+ch->size()+64);// Reserve space
        bit.reserve(bit.size()+ch->size()+64);  // Reserve space
        size_t idx{0};                  // Index into channel bitsream
        const size_t nbits=ch->size();  // Number of channel bits
        while (true)                    // For.... almost ever.....
        {                               // Process channel bits
          // ~~~~~~~~~~~~~~~~~~~~~~~~~~ //
          // Determine how many parity bits are available to this input step
          // ~~~~~~~~~~~~~~~~~~~~~~~~~~ //
          uint8_t hp1{0},hp2{0};        // Received parity bits
          uint8_t r1{0},r2{0};          // Received flags
          if (!ccfg.p34)                // Rate 1/2?
          {                             // Yes, so get both channel bits per input bit
            if (idx+1>=nbits)           // Not enough bits left?
              break;                    // Done processing
            r1=(*ch)[idx++]&1u;         // Get parity bit 1
            r2=static_cast<uint8_t>(((*ch)[idx++]&1u)^1u); // Get parity bit 2;
            hp1=1;                      // Mark parity bit 1 received
            hp2=1;                      // Mark parity bit 2 received
            ////Deb(" DecodeBits: Rate 1/2, r1=%d r2=%d",r1,r2);
          }
          else                          // Else it is 3/4 punctured
          {                             // 
            // ~~~~~~~~~~~~~~~~~~~~~~~~ //
            // Rate 3/4 puncturing phase: (2,1,1) bits per 3 input bits
            // ~~~~~~~~~~~~~~~~~~~~~~~~ //
            if (ph==0)                  // Is the phase counter 0?
            {                           // Yes, so get both parity bits
              if (idx+1>=nbits)         // Not engough bits left?
                break;                  // Done processing
              r1=(*ch)[idx++]&1u;       // Get parity bit 1
              r2=(*ch)[idx++]&1u;       // Get parity bit 2
              hp1=1;                    // Mark parity bit 1 received
              hp2=1;                    // Mark parity bit 2 received
              ////Deb(" DecodeBits: Rate 3/4, ph=%d r1=%d r2=%d",ph,r1,r2);
            }                           // Done with phase 0
            else if (ph==1)             // At phase 1 (t=2)?
            {                           // Yes, get only parity bit 1
              if (idx>=nbits)           // Not enough bits left?
                break;                  // Done processing
              r2=(*ch)[idx++]&1u;       // Get parity bit 1
              hp1=0;                    // Mark parity bit 1 received
              hp2=1;                    // Parity bit 2 not received
              ////Deb(" DecodeBits: Rate 3/4, ph=%d r1=%d",ph,r1);
            }                           // Done with phase 1
            else                        // Else phase 2
            {                           // Get only parity bit 2
              if (idx>=nbits)           // Not enough bits left?
                break;                  // Done processing
              r1=(*ch)[idx++]&1u;       // Get parity bit 2
              hp1=1;                    // Parity bit 1 not received
              hp2=0;                    // Mark parity bit 2 received
              //Deb(" DecodeBits: Rate 3/4, ph=%d r2=%d",ph,r2);
            }                           // Done with phase 2
          }                             // Done with 3/4 punctured code
          // ~~~~~~~~~~~~~~~~~~~~~~~~~~ //
          // One Viterbi step or phase recovers or corrects (one information bit)
          // ~~~~~~~~~~~~~~~~~~~~~~~~~~ //
          ViterbiStep(r1,r2,hp1,hp2);   // Perform one Viterbi time step
          // ~~~~~~~~~~~~~~~~~~~~~~~~~~ //
          // Advance puncturing phase
          // ~~~~~~~~~~~~~~~~~~~~~~~~~~ //
          if (ccfg.p34)                 // Rate 3/4?
            ph=(ph+1)%3;                // Advance puncturing phase
        }                               // Done processing channel bits
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        // Traceback to find best path
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        std::vector<uint8_t> urev{};    // Reversed input bits
        Traceback(&urev);               // Perform traceback
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        // Output in forward order
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        o->assign(urev.rbegin(),urev.rend());// Assign in forward order
        //Log("[END] Viterbi DecodeBits: output input bits=%zu",o->size());
      }                                 // ~~~~~~~~~~ DecodeBits ~~~~~~~~~~ //
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      // DecodeBitsWeighted — Soft-Decision (Soft-Input) Viterbi, and the road to SOVA.
      //
      // THEORY
      //   Classic Viterbi is Maximum-Likelihood Sequence Estimation (MLSE): it finds
      //   the single trellis path  u_hat = argmax_u P(r|u), minimizing the FRAME
      //   (sequence) error probability. It uses hard bits + a Hamming branch metric
      //   and emits a hard survivor with NO per-bit confidence.
      //
      //   This routine is the SOFT-INPUT variant. Each received bit carries a weight
      //   w in [0,1] (channel reliability) and the branch metric is a weighted
      //   mismatch COST (lower = better):
      //        BM(s->ns) = sum_i  w_i * [ c_i(s,ns) != r_i ]
      //   an unreliable bit (w->0) barely moves the path metric. This beats hard
      //   Viterbi (it exploits the analog margin) but the criterion is UNCHANGED:
      //   still MLSE, still a HARD output.
      //
      //   SOVA (Soft-Output Viterbi Algorithm, Hagenauer & Hoeher 1989) bolts a
      //   per-bit RELIABILITY onto the MLSE survivor, making Viterbi Soft-In/Soft-Out
      //   (SISO). At every ACS merge the discarded competitor has metric M_comp and
      //   the survivor M_surv; their margin
      //        Delta = M_comp - M_surv >= 0
      //   says how nearly the runner-up won. The bit APP magnitude is
      //        |L(u_k)| ~= min{ Delta_j : competitor at node j disagrees with the
      //                                   survivor on bit k },   sign(L)=hard decision.
      //   A small Delta at a node where the two paths disagree on bit k => bit k is
      //   unreliable.
      //
      //   RELATION TO BCJR (MAP): BCJR marginalizes over ALL paths (forward-backward)
      //   for the EXACT per-bit APP LLR and minimizes BIT error probability. SOVA
      //   keeps only the single best competitor, so it is a max-log-MAP-style
      //   approximation:  |L_SOVA| <= |L_BCJR|  (~0.5 dB weaker in iterative use),
      //   at ~2x hard-Viterbi cost vs BCJR's ~3-4x plus whole-block memory.
      //
      //     |                | criterion      | output      | cost        |
      //     | Viterbi (hard) | MLSE (seq)     | hard        | 1x          |
      //     | this (soft-in) | MLSE (seq)     | hard        | ~1x         |
      //     | SOVA           | MLSE + reliab. | hard + LLR  | ~2x         |
      //     | BCJR (MAP)     | MAP (per-bit)  | exact LLR   | ~3-4x + mem |
      //
      //   STATUS: this is now full SOVA. ViterbiStepWeighted captures the ACS margin
      //   Delta and the discarded competitor at every node; TracebackSova runs the
      //   reliability-update pass. Both classic SOVA outputs are exposed since they
      //   are computed together:
      //     rel -> per-bit reliability MAGNITUDE |L| (>=0), units of cost/BM_SCALE;
      //     llr -> signed classic-SOVA soft output L = sign(u_k)*|L|, L>0 favors
      //            bit 0 / L<0 favors bit 1 (same convention as BCJRDecoder, so it
      //            drops straight into the RS/frame-confidence path).
      //   + per-bit |L|, units of cost/BM_SCALE); pass rel == nullptr (default) for
      //   the classic soft-IN / hard-OUT MLSE traceback. Path metrics ARE
      //   renormalized every step (the minimum survivor is slid to zero), so
      //   frames of any length are safe. If you need the EXACT per-bit APP LLR, BCJRDecoder remains the
      //   better choice (SOVA's |L| is a max-log-MAP approximation, |L_SOVA| <= |L_BCJR|).
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //      
      inline void DecodeBitsWeighted (
       const std::vector<uint8_t>* ch,  // Input channel bits (hard 0/1)
        std::vector<float>* w,          // Weights for input channel bits
        std::vector<uint8_t>* const o,  // Output decoded input bits
        std::vector<float>* rel=nullptr,// Optional per-bit SOVA reliability MAGNITUDE |L| (>=0)
        std::vector<float>* llr=nullptr)// Optional signed classic-SOVA soft output L (L>0 favors bit 0)
      {                                 // ~~~~~~~~~~ DecodeBitsWeighted ~~~~~~~~~~ //
        if (ch==nullptr||w==nullptr||o==nullptr||ch->empty())
          return;
       if (ch->size()!=w->size())      // Size mismatch?
          return;                      // Nothing to do, size must match.
        o->clear();                    // Clear output buffer
        //Log("[BEGIN] Viterbi DecodeBitsWeighted: input channel bits=%zu",ch->size());
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        // Prepare survivor path storage
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        prev.reserve(prev.size()+ch->size()+64);// Reserve space
        bit.reserve(bit.size()+ch->size()+64);  // Reserve space
        size_t idx{0};                  // Index into channel bitsream
        const size_t nbits=ch->size();  // Number of channel bits
        while (true)                    // For.... almost ever.....
        {
          uint8_t r1{0},r2{0};          // Received parity bits
          uint8_t hp1{0},hp2{0};        // Received flags
          float w1{0.f},w2{0.f};        // Weights for received bits
          if (!ccfg.p34)                // Rate 1/2?
          {
            if (idx+1>=nbits)           // Not enough bits left?
              break;                    // Done processing
            r1=(*ch)[idx]&1u;            // Get parity bit 1
            w1=std::clamp((*w)[idx],0.f,1.f);// Get weight for parity bit 1
            idx++;                      // Advance index
            r2=static_cast<uint8_t>(((*ch)[idx]&1u)^1u); // Get parity bit 2
            w2=std::clamp((*w)[idx],0.f,1.f);// Get weight for parity bit 2
            idx++;                      // Advance index
            hp1=1;                      // Mark parity bit 1 received
            hp2=1;                      // Mark parity bit 2 received
            //Deb(" DecodeBitsWeighted: Rate 1/2, r1=%d w1=%.2f r2=%d w2=%.2f",r1,w1,r2,w2);
          }                             // Done with 1/2
          else
          {
            // ~~~~~~~~~~~~~~~~~~~~~~~~ //
            // Rate 3/4 puncturing phase: (0/1/2) bits per 3 input bits
            // ~~~~~~~~~~~~~~~~~~~~~~~~ //
            if (ph==0)                  // Phase counter zero?
            {
              if (idx+1>=nbits)         // Not engough bits left?
                break;                  // Done processing
              r1=(*ch)[idx]&1u;          // Get parity bit 1
              w1=std::clamp((*w)[idx],0.f,1.f);// Get weight for parity bit 1
              idx++;                    // Advance index
              r2=(*ch)[idx]&1u;          // Get parity bit 2
              w2=std::clamp((*w)[idx],0.f,1.f);// Get weight for parity bit 2
              idx++;                    // Advance index
              hp1=1;                    // Mark parity bit 1 received
              hp2=1;                    // Mark parity bit 2 received
              //Deb(" DecodeBitsWeighted: Rate 3/4, ph=%d r1=%d w1=%.2f r2=%d w2=%.2f",ph,r1,w1,r2,w2);
            }                           // Done with phase 0
            else if (ph==1)             // Phase 1?
            {                           // Yes
              if (idx>=nbits)           // Not enough bits left?
                break;                  // Done processing
              r2=(*ch)[idx]&1u;          // Get parity bit 1
              w2=std::clamp((*w)[idx],0.f,1.f);// Get weight for parity bit 1
              idx++;                    // Advance index
              hp1=0;                    // Mark parity bit 1 received
              hp2=1;                    // Parity bit 2 not received
              //Deb(" DecodeBitsWeighted: Rate 3/4, ph=%d r1=%d w1=%.2f",ph,r1,w1);
            }                           // Done with phase 1
            else                        // Else phase 2
            {                           //
              if (idx>=nbits)           // Not enough bits left?
                break;                  // Done processing
              r1=(*ch)[idx]&1u;          // Get parity bit 2
              w1=std::clamp((*w)[idx],0.f,1.f);// Get weight for parity bit 2
              idx++;                    // Advance index
              hp1=1;                    // Parity bit 1 not received
              hp2=0;                    // Mark parity bit 2 received
              //Deb(" DecodeBitsWeighted: Rate 3/4, ph=%d r2=%d w2=%.2f",ph,r2,w2);
            }                           // Done with phase 2
          }                             // Done with 3/4 punctured code
          ViterbiStepWeighted(r1,r2,hp1,hp2,w1,w2);// Perform one Viterbi time step
          if (ccfg.p34)                 // Rate 3/4?
            ph=(ph+1)%3;                // Advance puncturing phase
        }                               // Done processing channel bits
        if (rel!=nullptr||llr!=nullptr) // Soft output requested (magnitude and/or signed LLR)?
          TracebackSova(o,rel,30,llr);  // SOVA: hard bits + reliability |L| and/or signed LLR, forward order
        else                            // Hard output only
        {                               // Classic MLSE traceback
          std::vector<uint8_t> urev{};  // Reversed input bits
          Traceback(&urev);             // Perform traceback
          o->assign(urev.rbegin(),urev.rend());// Assign in forward order
        }                               // Done producing output
        //Log("[END] Viterbi DecodeBitsWeighted: output input bits=%zu",o->size());
      }                                 // ~~~~~~~~~~ DecodeBitsWeighted ~~~~~~~~~~ //
    public:
      // ~~~~~~~~~~~~~~~~~~~~~~~~~ //
      // Trellis structure for Viterbi Decoder
      // ~~~~~~~~~~~~~~~~~~~~~~~~~ //
      struct TrellisNode
      {
        private:
        uint8_t state{0};         // Current state
        uint8_t input{0};         // Input bit (0,1)
        uint8_t next[64][2];      // Next state for each input bit (0,1)
        uint8_t out[64][2][2];       // Output parity bits for each input bit (0,1)
        public:
        TrellisNode (void)=default;
        ~TrellisNode (void)=default;
        inline void Clear (void)
        {                         // ~~~~~~~~~~ Clear ~~~~~~~~~~ //
          for (int i=0;i<64;i++)  // For each state
          {                       //
            next[i][0]=0;         // Clear next state for input 0
            next[i][1]=0;         // Clear next state for input 1
            for (int j=0;j<2;j++)   // For each input bit
            {                     //
              out[i][j][0]=0;     // Clear output parity bit 1
              out[i][j][1]=0;     // Clear output parity bit 2
            }                     // Done for each input bit
          }                       //
        }                         // ~~~~~~~~~~ Clear ~~~~~~~~~~ //
        inline void Set (
          int state,              // Current state
          int input,              // Input bit (0,1)
          int nextstate,          // Next state
          int outbits)            // Output parity bits
        {                         // ~~~~~~~~~~ Set ~~~~~~~~~~ //
          next[state][input]=static_cast<uint8_t>(nextstate);// Set next state
          out[state][input][0]=static_cast<uint8_t>(outbits);  // Set output parity bits
        }                         // ~~~~~~~~~~ Set ~~~~~~~~~~ //
        inline void Get (
          int state,              // Current state
          int input,              // Input bit (0,1)
          int& nextstate,         // Next state (output)
          int& outbits) const     // Output parity bits (output)
        {                         // ~~~~~~~~~~ Get ~~~~~~~~~~ //
          nextstate=next[state][input];// Get next state
          outbits=out[state][input][0];  // Get output parity bits
        }                         // ~~~~~~~~~~ Get ~~~~~~~~~~ //
        inline uint8_t* GetNext (void)
        {
          return &this->next[0][0];
        }
        inline uint8_t* GetOut (void)
        {
          return &this->out[0][0][0];
        }
        inline uint8_t* GetNext (int state, int input)
        {
          return &this->next[state][input];
        }
        inline uint8_t* GetOut (int state, int input, int bit)
        {
          return &this->out[state][input][bit];
        }
        inline uint8_t GetCurrentState (void) const
        {
          return this->state;
        }
        inline void SetCurrentState (int state)
        {
          this->state=static_cast<uint8_t>(state);
        }
        inline uint8_t GetInputBit (void) const
        {
          return this->input;
        }
        inline void SetInputBit (int input)
        {
          this->input=static_cast<uint8_t>(input);
        }
        inline uint8_t Parity (uint32_t bs) // Bit sequence to parity check
        {                               // ~~~~~~~~~~ Parity ~~~~~~~~~~ //                                
          bs^=bs>>16u;                  // Fold upper bits
          bs^=bs>>8u;                   // Fold again
          bs^=bs>>4u;                   // Fold again
          bs^=bs>>2u;                   // Fold again
          bs^=bs>>1u;                   // Fold again
          return static_cast<uint8_t>(bs&0x1u);// Return parity bit
        }                               // ~~~~~~~~~~ Parity ~~~~~~~~~~ //
        inline void BuildTrellis (void)
        {                               // ~~~~~~~~~~~~ BuildTrellis ~~~~~~~~~~~~ //
          Clear();                      // Clear existing trellis
          const uint32_t g1=0b1111001u; // Generator 1
          const uint32_t g2=0b1011011u; // Generator 2
          for (uint8_t s=0;s<64;++s)    // For each state
          {                             // and for each input bit
            for (uint8_t b=0;b<2;++b)   // and for each input bit.....
            {                           // Compute Trellis diagram
              const int zprm=((s<<1)|b)&0x7F; // Shift register with input bit
              next[s][b]=static_cast<uint8_t>(zprm&0x3F);// Next state (6 LSBs)
              out[s][b][0]=static_cast<uint8_t>(Parity(static_cast<uint32_t>(zprm)&g1));// Expected C1 (G1)
              // ~~~~~~~~~~~~~~~~~~~~~~ //
              // CCSDS 131.0-B-5 sec 3.3.1(5): G2 output inversion. Mirrors the encoder
              // so the trellis expects the inverted C2 for rate 1/2 AND 3/4
              // (the puncture metric reads this same table entry)
              // ~~~~~~~~~~~~~~~~~~~~~~ //
              out[s][b][1]=static_cast<uint8_t>(Parity(static_cast<uint32_t>(zprm)&g2));// RAW
              // //Deb(" BuildTrellis: State %02X Input %d -> Next %02X Output [%d,%d]",
              //    s,b,next[s][b],out[s][b][0],out[s][b][1]);
            }                           // Done for each input bit
          }                             // Done for each state
        }                               // ~~~~~~~~~~~~ BuildTrellis ~~~~~~~~~~~~ //
        inline void Print (void) const
        {                               // ~~~~~~~~~~ Print ~~~~~~~~~~ //
          for (uint8_t s=0;s<64;++s)    // For each state
          {                             //
            for (uint8_t b=0;b<2;++b)   // For each input bit
            {                           //
              uint8_t ns=next[s][b];    // Next state
              uint8_t ob1=out[s][b][0]; // Output bit 1
              uint8_t ob2=out[s][b][1]; // Output bit 2
              //std::printf(" State %02X Input %d -> Next %02X Output [%d,%d]\n",
              //  s,b,ns,ob1,ob2);
              Log(" State %02X Input %d -> Next %02X Output [%d,%d]",
                 s,b,ns,ob1,ob2);
            }                           // Done for each input bit
          }                             // Done for each state
        }                               // ~~~~~~~~~~ Print ~~~~~~~~~~ //        
      } trellis{};      
    private:
      ConvConfig ccfg{};                // Convolutional Codec conf parameters
      int ph{0};                        // Puncturing phase counter
      static constexpr uint16_t INF=0x3FFFu;// big and safe for short frames
      static constexpr uint16_t BM_SCALE=256u; // Branch metric scaling factor
      static constexpr float BM_SCALE_F32=static_cast<float>(BM_SCALE);
      uint16_t pm[64]{};                 // Path metrics
      size_t steps{0};                  // Number of steps Viterbi has gone through
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      // Survivor paths: one row per step, 64 colums (state)
      // Each row stores predecessor state and the input bit taken to enter that state
      // and lead us to that path at that step.
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      std::vector<std::array<uint8_t,64>> prev{};   // Survivor paths
      std::vector<std::array<uint8_t,64>> bit{};    // Input bits
      std::vector<std::array<uint16_t,64>> dmarg{}; // SOVA: ACS margin Delta per step/state
      std::vector<std::array<uint8_t,64>> prevc{};  // SOVA: competitor predecessor state
      std::vector<std::array<uint8_t,64>> bitc{};   // SOVA: competitor input bit
      // Log object

      // One ACS (Add path weights, Compare, Select best path) step
      inline void ViterbiStep (
        uint8_t r1,                     // Received parity bit 1
        uint8_t r2,                     // Received parity bit 2
        uint8_t isp1,                   // Do we have parity bit 1?
        uint8_t isp2)                   // Do we have parity bit 2?
      {                                 // ~~~~~~~~~~ ViterbiStep ~~~~~~~~~~ //
        //if (lg)
        //Deb(" [BEGIN] ViterbiStep: r1=%d r2=%d isp1=%d isp2=%d",r1,r2,isp1,isp2);
        // Get new path metrics
        uint16_t qm[64]{};              // New path metrics
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        // Store survivors for this step
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        uint8_t ps[64];                 // Predecessor state (prev row)
        uint8_t ib[64];                 // Input bit that led to this survivor state.
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        // Init new path metrics to Inf
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        for (int s=0;s<64;++s)          // For each state
          qm[s]=INF;                    // Init to Inf
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        // For each state, compute branch metrics for input bits 0 and 1
        // and update path metrics accordingly
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        for (int s=0;s<64;++s)          // For each state
        {                               // Viterbi Add path weigths, compare and select.....
          const uint64_t base=static_cast<uint64_t>(pm[s]);    // Base path metric for this state
          if (base>=INF)                // Invalid path?
            continue;                   // Skip it
          // ~~~~~~~~~~~~~~~~~~~~~~~~~~ //
          // Try two possibilities: input bit 0 and 1 and their respective outputs
          // ~~~~~~~~~~~~~~~~~~~~~~~~~~ //
          for (int b=0;b<2;++b)         // For each input bit
          {                             //
            const int ns=static_cast<int>(trellis.GetNext(s,b)[0]); // Next state
            const uint8_t* e1=trellis.GetOut(s,b,0); // Output bits for this input
            const uint8_t* e2=trellis.GetOut(s,b,1); // Output bits for this input
            // ~~~~~~~~~~~~~~~~~~~~~~~~ //
            // Compute Hamming distance (branch metric) from received bits to
            // expected output bits for this transition
            // ~~~~~~~~~~~~~~~~~~~~~~~~ //
            uint8_t bm{0};              // Branch metric
            if (isp1!=0)                // Parity bit 1 present?
              bm+=(*e1!=r1);            // Done computing branch metric
            if (isp2!=0)                // Parity bit 2 present?
              bm+=(*e2!=r2);            // Done computing branch metric
            //if (lg)
            //  //Deb(" ViterbiStep: State %02X Input %d -> Next %02X Expected [%d,%d] Received [%d,%d] BM=%d",
            //    s,b,ns,*e1,*e2,
            //    (isp1!=0)? r1:0xFF,
            //    (isp2!=0)? r2:0xFF,
             //   bm);
            // ~~~~~~~~~~~~~~~~~~~~~~~~ //
            // Compute new path metric for this transition
            // ~~~~~~~~~~~~~~~~~~~~~~~~ //
            const uint16_t npm=static_cast<uint16_t>(base+bm);// New path metric
            // ~~~~~~~~~~~~~~~~~~~~~~~~ //
            // Compare with existing path metric for next state
            // ~~~~~~~~~~~~~~~~~~~~~~~~ //
            if (npm<qm[ns])             // Better path?
            {                           // Yes, so update optimal path
              qm[ns]=npm;               // Update path metric
              ps[ns]=static_cast<uint8_t>(s);// Store predecessor state
              ib[ns]=static_cast<uint8_t>(b);// Store input bit that led to this state
              //if (lg)
              //  //Deb(" ViterbiStep: Update State %02X New PM=%d Prev State %02X Input Bit %d",
              //    ns,npm,s,b);
            }                           // Done updating optimal path
          }                             // Done for each input bit
        }                               // Done for each state
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        // Comment and renormalize.
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        uint16_t mn=INF;                // Running minimum over valid survivors
        for (int s=0;s<64;++s)          // For each state
        {                               // Commit and track minimum
         pm[s]=qm[s];                   // Commit new path metric
         if (pm[s]<mn)                  // New minimum
           mn=pm[s];                    // Yes so save it
        }                               // Done commit
        if (mn>0u&&mn<INF)              // Anything to subtract and valid path metric?
        {                               // Yes
         for (int s=0;s<64;++s)         // For each state
           if (pm[s]<INF)               // Valid path (not the sentinnel)?
             pm[s]=static_cast<pm[s]-mn); // Yes, slide the whole metric floor.
        }                               // Done sliding metric floor.
        prev.emplace_back();            // Add new row for predecessor states
        bit.emplace_back();             // Add new row for input bits
        for (int s=0;s<64;++s)          // For each state
        {                               // Store last step survivors and input bits
          prev.back()[s]=ps[s];         // Predecessor state
          bit.back()[s]=ib[s];          // Input bit that led to this state
        }                               // Done for each state
        ++steps;                        // Advance number of steps
        //if (lg)
        //  //Deb(" [END] ViterbiStep: Completed step %zu",steps);
      }                                 // ~~~~~~~~~~ ViterbiStep ~~~~~~~~~~ //      
      inline void ViterbiStepWeighted (
        uint8_t r1,                     // Received parity bit 1
        uint8_t r2,                     // Received parity bit 2
        uint8_t isp1,                   // Do we have parity bit 1?
        uint8_t isp2,                   // Do we have parity bit 2?
        float w1,                       // Weight for parity bit 1
        float w2)                       // Weight for parity bit 2
      {                                 // ~~~~~~~~~~ ViterbiStepWeighted ~~~~~~~~~~ //
        uint16_t qm[64];                // Best path metrics
        std::fill(std::begin(qm),std::end(qm),INF); // Init new path metrics to Inf
        uint16_t qm2[64];               // Best path metrics 2 for competitor (SOVA)
        std::fill(std::begin(qm2),std::end(qm2),INF); // Init new path metrics to Inf
        uint8_t ps[64]{};               // Survivor predecessor states
        uint8_t ib[64]{};               // Survivor input bits
        uint8_t psc[64]{};              // Competitor predecessor states (SOVA)
        uint8_t ibc[64]{};              // Competitor input bits (SOVA)
        const uint16_t w1s=static_cast<uint16_t>(std::clamp(w1,0.f,1.f)*BM_SCALE+0.5f); // Scaled weight for bit 1
        const uint16_t w2s=static_cast<uint16_t>(std::clamp(w2,0.f,1.f)*BM_SCALE+0.5f); // Scaled weight for bit 2
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        // For each state, compute branch metrics for input bits 0 and 1
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        for (int s=0;s<64;++s)          // For each state
        {                               // Viterbi add path wieghts, compare and select best path
          const uint16_t base=pm[s];    // Base path metric for this state
          if (base>=INF)                // Invalid path?
            continue;                   // Skip it
          for (int b=0;b<2;++b)         // For each input bit
          {                             // Walk trellis
            const int ns=static_cast<int>(trellis.GetNext(s,b)[0]); // Next state
            const uint8_t* e1=trellis.GetOut(s,b,0); // Expected output bits for this input
            const uint8_t* e2=trellis.GetOut(s,b,1); // Expected output bits for this input
            uint16_t bm{0};             // Branch metric
            if (isp1!=0)                // Parity bit 1 present?
              bm+=((*e1!=r1)?w1s:0);    // Add scaled weight if bit mismatch
            if (isp2!=0)                // Parity bit 2 present?
              bm+=((*e2!=r2)?w2s:0);    // Add scaled weight if bit mismatch
            const uint16_t npm=static_cast<uint16_t>(std::min<uint32_t>(INF,base+bm)); // New path metric
            if (npm<qm[ns])             // Better path metric for next state?
            {                           // Yes: new survivor, old survivor demotes to competitor
              qm2[ns]=qm[ns];           // Old best becomes the competitor metric
              psc[ns]=ps[ns];           // Old survivor predecessor becomes competitor predecessor
              ibc[ns]=ib[ns];           // Old survivor input bit becomes competitor input bit
              qm[ns]=npm;               // Update best (survivor) path metric
              ps[ns]=static_cast<uint8_t>(s); // Survivor predecessor STATE (not a metric)
              ib[ns]=static_cast<uint8_t>(b); // Update survivor input bit
            }                           // Done updating survivor
            else if (npm<qm2[ns])       // Else bettern than competitor?
            {                           // Update new competitor
              qm2[ns]=npm;              // Update competitor path metric
              psc[ns]=static_cast<uint8_t>(s); // Update competitor predecessor state
              ibc[ns]=static_cast<uint8_t>(b); // Update competitor input bit
            }                           // Done updating competitor
            else                        // NO OP IT
            {                           // Competitor not better than survivor, do nothing
              // NO-OP
            }                           // Done comparing with survivor and competitor
          }                             // Done for each input bit
        }                               // Done for each state
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        // Commit new path metrics and push survivors for traceback
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        uint16_t mn=INF;                // Running minumum over valid survivors
        for (int s=0;s<64;++s)          // For each state
        {                               // Commit and track min
          pm[s]=qm[s];                  // Commit new path metric with best path metric
          if (pm[s]<qm[s])              // New min?
            mn=pm[s];                   // Yes, store that.
        }                               // Done committing
        if (mn>0u&&mn<INF)              // Anything to subtract and valid path metric?
        {                               // Yes
          for (int s=0;s<64;++s)        // For each state
            if (pm[s]<INF)              // Valid path?
              pm[s]=static_cast<uint16_t>(pm[s]-mn); // Slide the whole metric floor
        }
        prev.emplace_back();            // Add new row for predecessor states
        bit.emplace_back();             // Add new row for input bits
        dmarg.emplace_back();           // Add new row for SOVA margins
        prevc.emplace_back();           // Add new row for SOVA competitor predecessor states
        bitc.emplace_back();            // Add new row for SOVA competitor input bits
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        // Update survivor and competitor storage for traceback
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        for (int s=0;s<64;++s)          // For each state
        {                               // Store survivor and competitor info for this step/state
          prev.back()[s]=ps[s];         // Survivor predecessor state
          bit.back()[s]=ib[s];          // Survivor input bit that led to THIS state
          // ~~~~~~~~~~~~~~~~~~~~~~~~~~ //
          // Compute SOVA margin delta = competitor - survivor (>=0). If either
          // path is invalid (no competitor merged), the margin is INF (certain).
          // ~~~~~~~~~~~~~~~~~~~~~~~~~~ //
          const uint16_t d=(qm[s]>=INF||qm2[s]>=INF)?INF:static_cast<uint16_t>(qm2[s]-qm[s]); // Competitor margin
          dmarg.back()[s]=d;            // Store SOVA margin for this step/state
          prevc.back()[s]=psc[s];       // Store competitor predecessor state
          bitc.back()[s]=ibc[s];        // Store competitor input bit
        }                               // Done for each state
        ++steps;                        // Advance number of steps
      }                                 // ~~~~~~~~~~ ViterbiStepWeighted ~~~~~~~~~~ //
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      // Traceback to find best path through trellis and output decoded bits in reverse order
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      inline void Traceback (
        std::vector<uint8_t>* const urev) // Reversed output bits
      {                                 // ~~~~~~~~~~ Traceback ~~~~~~~~~~~~ //
        if (urev==nullptr||steps==0)              // Bad arg?
          return;                       // Nothing to do
       // if (lg)
       //   //Log("[BEGIN] Viterbi Traceback");
        urev->clear();                  // Clear output buffer
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        // Find state with best path metric
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        uint16_t bstm=INF;              // Best path metric
        int bsts{0};                    // State with best path metric
        for (int s=0;s<64;++s)          // For each state
        {                               // Find best path metric
          if (pm[s]<bstm)               // Do we have a better path metric?
          {                             // Yes
            bstm=pm[s];                 // Update best path metric
            bsts=s;                     // Update best state
          }                             // Done checking this state
        }                               // Done for all states
        //Log(" Viterbi Traceback: Best State=%02X Best Path Metric=%d",bsts,bstm);
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        // Traceback through survivor paths
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        const int TB{30};               // Traceback depth
        int tbstps=(steps>TB)?(steps-TB):0;// Traceback start step
        int ste=bsts;                   // Best state to start traceback
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        // Perform a walkback through all steps, but emit all decisions
        // only after reaching the traceback depth.
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        for (int t=static_cast<int>(steps)-1;t>=0;--t)// For each step in diagram....
        {                               // Trace backwards....
          const uint8_t u=bit[static_cast<size_t>(t)][ste];// Input bit that led to this state
          const int ps=prev[static_cast<size_t>(t)][ste];// Predecessor state
          if (t<tbstps)                 // Beyond traceback depth?
            urev->push_back(u);         // Yes, so keep reliable tail first
          else                          // Else not yet at traceback depth
            urev->push_back(u);         // Keep all bits 
          ste=ps;                       // Move to predecessor state
         // if (lg)
          //  //Deb(" [END] Viterbi Traceback: Step %d Current State %02X Input Bit %d Predecessor State %02X",
          //    t,ste,u,ps);
        }                               // Done for all steps
      }                                 // ~~~~~~~~~~ Traceback ~~~~~~~~~~~~ //
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      // SOVA reliability update pass: after traceback, we have the survivor path and its decisions
      // for all steps. For each step, we can look at the competitor path that lost at that step and see if it disagreed with the survivor on the decoded bit. If so,
      // we can compute the metric difference Delta between the survivor and competitor at that step and use it to update the reliability of the decoded bit. The larger the Delta, the more reliable the decoded bit is (the more the survivor won by).
      // U controls how many steps to look back for reliability updates: U=0 means only the last bit, U=1 means the last two bits, and so on, up to U
      // =steps means all bits. In practice, a small U (e.g. 5-10) may be sufficient to capture most of the reliability information while keeping complexity manageable.
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      inline void TracebackSova (
        std::vector<uint8_t>* out,      // Output decoded bits (0/1)
        std::vector<float>* rel,        // Output per-bit reliability MAGNITUDE |L| (>=0); nullptr to skip
        int U=30,                       // Traceback depth for reliability update (0=only last bit, U=all bits)
        std::vector<float>* llr=nullptr)// Optional signed classic-SOVA soft output L (L>0 favors bit 0); nullptr to skip
      {                                 // ~~~~~~~~~~ TracebackSova ~~~~~~~~~~ //
        if (out==nullptr)               // Bad args?
          return;                       // Nothing to do
        out->clear();                   // Clear output bits buffer
        if (rel!=nullptr)               // We have reliability buffer?
          rel->clear();                 // Clear reliability buffer
        if (llr!=nullptr)               // We have signed soft-output buffer?
          llr->clear();                 // Clear signed LLR buffer
        if (steps==0) return;           // No steps, nothing to do
        uint16_t bstm=INF;              // Best path metric
        int bsts=0;                     // Final state with best path metric
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        // Find the best path metric
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        for (int s=0;s<64;++s)          // For each state
        {                               // Find best path metric
          if (pm[s]<bstm)               // Do we have a better path metric?
          {                             // Yes
            bstm=pm[s];                 // Update best path metric
            bsts=s;                     // Update best state
          }                             // Done checking this state
        }                               // Done for all states
        const int tst=static_cast<int>(steps);// Total number of steps
        std::vector<int> mls(tst);      // Maximum-Likelihood state per time (fwd)
        std::vector<uint8_t> mlb(tst);  // Maximum-Likelihood bit per time 
        {                               // Lock scope scope for ML path computation
          int s=bsts;                   // Start from best state at end
          for (int t=tst-1;t>=0;--t)    // For each step in the diagram, backwards
          {                             // Trace back survivor path to find ML decisions and states
            mls[static_cast<size_t>(t)]=s; // Store ML state at this time
            mlb[static_cast<size_t>(t)]=bit[static_cast<size_t>(t)][s]; // Store ML bit at this time
            s=prev[static_cast<size_t>(t)][s]; // Move to predecessor state
          }                             // Done tracing back survivor path
        }                               // Done with ML path computation
        std::vector<uint16_t> ell(tst,INF); // Relaibility L (cost metric units)
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        // Compute reliability for each decision
        // along the survivor path by looking
        // at the competitor path that lost to the survivor
        // at each time step.
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        for (int t=tst-1;t>=0;--t)      // For each step in the diagram, backwards
        {                               // Compute reliability for ML decision at THIS time
          size_t idx=static_cast<size_t>(t); // Index for current time step
          const int s=mls[idx];         // ML state at this time
          const uint16_t d=dmarg[idx][s]; // Metric margin for ML path at this time
          if (d>=INF)                   // Invalid margin?
            continue;                   // Skip reliability update for this time
          if (bitc[idx][s]!=mlb[idx])   // Competitor bit disagrees with ML bit?
            ell[idx]=std::min<uint16_t>(ell[idx],d); // Update reliability metric for this time
          int cs=prevc[idx][s];         // Competitor state at this time
          // ~~~~~~~~~~~~~~~~~~~~~~~~~~ //
          // Trace back competitor path for U steps
          // or until path re-merges with survivor path,
          // whichever comes first, and update reliability for any disagreements
          // along the way.
          // ~~~~~~~~~~~~~~~~~~~~~~~~~~ //
          for (int i=t-1;i>=0&&i>=t-U;--i) // For previous U steps
          {                             // Walk competitor back through trellisto the re-merge
            if (cs==mls[i])             // Paths have re-merged by this time?
              break;                    // Done with this competitor path
            if (bit[i][cs]!=mlb[i])     // Competitor bit disagrees with ML bit?
              ell[idx]=std::min<uint16_t>(ell[idx],d); // Update reliability metric for this time
            cs=prev[static_cast<size_t>(i)][cs]; // Move competitor to predecessor state
          }                             // Done walking to prev U states
        }                               // Done with this time step
        out->resize(tst);               // Resize output bits buffer
        if (rel!=nullptr)               // We have reliability buffer?
          rel->resize(tst);             // Resize reliability buffer
        if (llr!=nullptr)               // We have signed soft-output buffer?
          llr->resize(tst);             // Resize signed LLR buffer
        for (int t=0;t<tst;++t)         // For each time step
        {                               // Recollect maximally likely determined bits
          size_t idx=static_cast<size_t>(t); // Index for this time step
          (*out)[idx]=mlb[idx];         // Output ML bit for this time step
          const float mag=(ell[idx])/BM_SCALE_F32; // Reliability MAGNITUDE |L| (cost/BM_SCALE)
          if (rel!=nullptr)             // We have a reliability buffer?
            (*rel)[idx]=mag;            // Output reliability magnitude for this time step
          if (llr!=nullptr)             // We have a signed soft-output buffer?
            (*llr)[idx]=(mlb[idx]?-mag:mag); // Classic SOVA LLR: +|L| for bit 0, -|L| for bit 1
        }                               // Done for each time step
      }                                 // ~~~~~~~~~~ TracebackSova ~~~~~~~~~~ //
  };
}
