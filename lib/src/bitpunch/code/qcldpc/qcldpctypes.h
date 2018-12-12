/*
This file is part of BitPunch
Copyright (C) 2015 Frantisek Uhrecky <frantisek.uhrecky[what here]gmail.com>
Copyright (C) 2015 Andrej Gulyas <andrej.guly[what here]gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef QCLDPCTYPES_H
#define QCLDPCTYPES_H
#include <bitpunch/math/gf2types.h>
#include <bitpunch/debugio.h>
#include <stdlib.h>

/**
 * QC-LDPC McEliece code matrices
 */
typedef struct _BPU_T_Qcldpc_Spec {
	BPU_T_GF2_QC_Matrix G_masked; ///< masked Generating Matrix
	BPU_T_GF2_QC_Matrix S; ///< dense masking matrix S
  	BPU_T_GF2_Sparse_Qc_Matrix H; ///< parity-check matrix
	BPU_T_GF2_Sparse_Qc_Matrix Q_sparse; ///< sparse masking matrix Q
	uint16_t m; ///< length of a generating polynomial
	uint16_t n0; ///< number of cyclic matrices H
	uint16_t w; ///< weight of parity-check matrix row	
}BPU_T_Qcldpc_Spec;

/**
  QC-LDPC McEliece input code params
  */
typedef struct _BPU_T_Qcldpc_Params {
	uint16_t m; ///< length of a generating polynomial
	uint16_t n0; ///< number of cyclic matrices H
	uint16_t w; ///< weight of parity-check matrix row
	uint16_t t; ///< weight of error vector	
}BPU_T_Qcldpc_Params;

/**
 * Free QC-LDPC McEliece code matrices
 * @param spec pointer to structure
 */
void BPU_qcldpcFreeSpec(BPU_T_Qcldpc_Spec *spec);

/**
 * Allocate memory for QC-LDPC code params. After work you have to free memory using call BPU_qcldpcFreeParams
 * @param  params pointer to structure
 * @param  m      length of a generating polynomial
 * @param  n0     number of cyclic matrices H
 * @param  w      weight of parity-check matrix row
 * @param  t      weight of error vector
 * @return        0 if OK, else error
 */
int BPU_qcldpcInitParams(BPU_T_Qcldpc_Params **params, const uint16_t m, const uint16_t n0, const uint16_t w, const uint16_t t);

/**
 * Free memory for QC-LDPC code params.
 * @param params pointer to structure
 */
void BPU_qcldpcFreeParams(BPU_T_Qcldpc_Params **params);

#endif // QCLDPCTYPES_H
