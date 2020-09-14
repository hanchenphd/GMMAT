#ifndef READBGEN_H
#define READBGEN_H

#include <R.h>
#include <Rcpp.h>
using namespace Rcpp;
#include <fstream>
#include <cstring>
#include <cstdio>
#include <vector>

typedef unsigned int uint;
typedef unsigned char uchar;
typedef unsigned short ushort;

using namespace std;
using namespace Rcpp;


void Bgen13GetTwoVals(const unsigned char* prob_start, uint32_t bit_precision, uintptr_t offset, uintptr_t* first_val_ptr, uintptr_t* second_val_ptr) {
  
  switch(bit_precision) {
  case 8:
    *first_val_ptr  = prob_start[0];
    prob_start += offset;
    *second_val_ptr = prob_start[0];
    break;
  case 16:
    *first_val_ptr  = prob_start[0]|(prob_start[1]<<8);
    prob_start += offset;
    *second_val_ptr = prob_start[0]|(prob_start[1]<<8);
    break;
  case 24:
    *first_val_ptr  = prob_start[0]|(prob_start[1]<<8)|(prob_start[2]<<16);
    prob_start += offset;
    *second_val_ptr = prob_start[0]|(prob_start[1]<<8)|(prob_start[2]<<16);
    break;
  case 32:
    *first_val_ptr  = prob_start[0]|(prob_start[1]<<8)|(prob_start[2]<<16)|(prob_start[3]<<24);
    prob_start += offset;
    *second_val_ptr = prob_start[0]|(prob_start[1]<<8)|(prob_start[2]<<16)|(prob_start[3]<<24);
    break;
  }
  
}



extern "C" 
{
  SEXP bgenHeader(SEXP bgenfile_in){
    string bgenfile = Rcpp::as<string>(bgenfile_in);
    
    FILE* bStream;
    bStream = fopen(bgenfile.c_str(), "rb");
    if (!bStream) { 
      Rcout << "Error reading BGEN file: " << bgenfile << "\n"; return R_NilValue; 
    }
    
    uint offset;   if(!fread(&offset, 4, 1, bStream)) { Rcout << "Error reading BGEN file: Cannot read offset value in header block. \n"; return R_NilValue;}
    uint L_H;      if(!fread(&L_H,    4, 1, bStream)) { Rcout << "Error reading BGEN file: Cannot read length of header block. \n"; return R_NilValue;}
    uint Mbgen;    if(!fread(&Mbgen,  4, 1, bStream)) { Rcout << "Error reading BGEN file: Cannot read number of variants in header block. \n"; return R_NilValue;}
    uint Nbgen;    if(!fread(&Nbgen,  4, 1, bStream)) { Rcout << "Error reading BGEN file: Cannot read number of samples in header block. \n"; return R_NilValue;}
    char magic[4]; if(!fread(magic,   1, 4, bStream)) { Rcout << "Error reading BGEN file: Cannot read number magic bytes in header block. \n"; return R_NilValue;}
    fseek(bStream, L_H - 20, SEEK_CUR);
    uint flags;    if(!fread(&flags,  4, 1, bStream)) { Rcout << "Error reading BGEN file: Cannot read flag in header block. \n"; return R_NilValue;}
    
    if (!(magic[0] == 'b' && magic[1] == 'g' && magic[2] == 'e' && magic[3] == 'n')){
      Rcout << "Error reading BGEN file: BGEN file's magic number does not match 'bgen'.\n"; return R_NilValue;
    }
    uint Compression = flags & 3;
    uint Layout      = (flags >> 2) & 0xf;
    uint SampleID    = flags >> 31;
    
    
    if (Compression > 2U){
      Rcout << "Error reading BGEN file: BGEN compression flag (" << Compression << ") should have a value of 0, 1 or 2.\n"; return R_NilValue;
    }
    if (Layout < 1U || Layout > 2U){
      Rcout << "Error reading BGEN file: BGEN layout flag (" << Layout<< ") should  have a value of 1 or 2.\n"; return R_NilValue;
    }
    if (SampleID > 1U){
      Rcout << "Error reading BGEN file: BGEN sample identifier flag (" << SampleID << ") should have a value of 0 or 1.\n"; return R_NilValue;
    }
    
    
    uint maxLA = 65536;
    vector<string> samVec;
    if (SampleID == 1){
      uint LS;  if(!fread(&LS,  4, 1, bStream)) {  Rcout << "Error reading BGEN file: Cannot read length of sample block value.\n"; return R_NilValue; }
      uint Nid; if(!fread(&Nid, 4, 1, bStream)) {  Rcout << "Error reading BGEN file: Cannot read the number of sample identifiers.\n"; return R_NilValue;}
      
      if (Nid != Nbgen) {
        Rcout << "Error reading BGEN file: Number of sample identifiers (" << Nid<< ") does not match number of samples specified in BGEN header block (" << Nbgen << ").\n"; return R_NilValue;
      }
      
      samVec.resize(Nid);
      char* samID = new char[maxLA + 1];
 
      for (uint n = 0; n < Nid; n++) {
          ushort LSID;
          if (!fread(&LSID, 2, 1, bStream)) { Rcout << "Error reading BGEN file: Cannot read in sample identifiers.\n"; return R_NilValue; }
          if (!fread(samID, 1, LSID, bStream)) { Rcout << "Error reading BGEN file: Cannot read in sample identifiers.\n"; return R_NilValue; }
          samID[LSID] = '\0';
          samVec[n] = string(samID);
      }
      delete[] samID;
    }
    

    fclose(bStream);
    return(Rcpp::List::create(Named("offset") = offset,
                              Named("M") = Mbgen,
                              Named("N") = Nbgen,
                              Named("CompressionFlag") = Compression,
                              Named("LayoutFlag") = Layout,
                              Named("SampleIdFlag") = SampleID,
                              Named("SampleIds") = samVec));
    
  }
  
  SEXP getVariantPos(SEXP bgenfile_in, SEXP offset_in, SEXP mbgen_in, SEXP nbgen_in, SEXP compression_in, SEXP layout_in, SEXP cores_in){
    
    string bgenfile  = Rcpp::as<string>(bgenfile_in);
    uint compression = Rcpp::as<uint>(compression_in);
    uint layout      = Rcpp::as<uint>(layout_in);
    uint mbgen       = Rcpp::as<uint>(mbgen_in);
    uint nbgen       = Rcpp::as<uint>(nbgen_in);
    uint offset      = Rcpp::as<uint>(offset_in);
    uint threads     = Rcpp::as<uint>(cores_in);

    uint maxLA = 65536;
    uint maxLB = 65536;
    char* snpID   = new char[maxLA + 1];
    char* rsID    = new char[maxLA + 1];
    char* chrStr  = new char[maxLA + 1];
    char* allele1 = new char[maxLA + 1];
    char* allele0 = new char[maxLB + 1]; 
    
    vector<uint> begin(threads);
    vector<uint> end(threads);
    vector<long long unsigned int> bgenVariantPos(threads);
    
      
    for (uint i = 0; i < threads; i++) {
      begin[i] = floor((mbgen / threads) * i);
      
      if ((i + 1) == (threads)) {
        end[i] = mbgen;
      }
      else {
        end[i] = floor(((mbgen / threads) * (i + 1)));
      }
    }
      
      
    FILE* fin = fopen(bgenfile.c_str(), "rb");
    fseek(fin, offset + 4, SEEK_SET);
      
    uint index = 0;
 

    for (uint snploop = 0; snploop < mbgen; snploop++) {
      int ret;
      if (snploop == begin[index]) {
        bgenVariantPos[index] = ftell(fin);
        index++;
        
        if (index == (begin.size())) {
            break;
        }
      }
        
      if (layout == 1) {
        uint Nid; ret = fread(&Nid, 4, 1, fin);  
        if (Nid != nbgen) {
          Rcout << "Error reading bgen file: Number of samples with genotype probabilties (" << Nid << ") does not match number of samples in BGEN file (" << nbgen << ").\n";
          return R_NilValue;
        }
      }
        
        
      ushort LS; ret = fread(&LS, 2, 1, fin);
      ret = fread(snpID, 1, LS, fin); snpID[LS] = '\0';
      
      ushort LR; ret = fread(&LR, 2, 1, fin);
      ret = fread(rsID, 1, LR, fin); rsID[LR] = '\0';
        
      ushort LC; ret = fread(&LC, 2, 1, fin);
      ret = fread(chrStr, 1, LC, fin); chrStr[LC] = '\0';
       
      uint physpos; ret = fread(&physpos, 4, 1, fin);
      
      ushort LKnum;
      if (layout == 2) {
        ret = fread(&LKnum, 2, 1, fin); 
        if (LKnum != 2) {
          Rcout << "Error reading bgen file: " << string(snpID) << " is a non-bi-allelic variant with " << LKnum << " alleles. Please filter these variants for now.";
          return R_NilValue;
        }
      }
      
      uint LA; ret = fread(&LA, 4, 1, fin);
      ret = fread(allele1, 1, LA, fin); allele1[LA] = '\0';
      
      uint LB; ret = fread(&LB, 4, 1, fin);
      ret = fread(allele0, 1, LB, fin); allele0[LB] = '\0';
      
      if (layout == 2) {
        if (compression > 0) {
          uint zLen; ret = fread(&zLen, 4, 1, fin);
          fseek(fin, 4 + zLen - 4, SEEK_CUR);
        }
        else {
          uint zLen; ret = fread(&zLen, 4, 1, fin);
          fseek(fin, zLen, SEEK_CUR);
        }
      }
      else {
        if (compression == 1) {
          uint zLen; ret = fread(&zLen, 4, 1, fin);
          fseek(fin, zLen, SEEK_CUR);
        }
        else {
          fseek(fin, 6 * nbgen, SEEK_CUR);
        }
        
      }

      if (ret == 0) { Rcout << "Error reading BGEN file : Cannot parse variant block.\n"; return R_NilValue; }
    }
    
    delete[] snpID;
    delete[] rsID;
    delete[] chrStr;
    delete[] allele0;
    delete[] allele1;

    fclose(fin);
    return(Rcpp::List::create(Named("begin") = begin,
                              Named("end") = end,
                              Named("pos") = bgenVariantPos));
  }
}

#endif



