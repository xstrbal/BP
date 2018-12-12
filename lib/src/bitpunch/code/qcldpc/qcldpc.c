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
#include <bitpunch/bitpunch.h>
#include "qcldpc.h"
#include <bitpunch/errorcodes.h>

#ifdef BPU_CONF_ENCRYPTION
int BPU_mecsQcldpcEncode1(BPU_T_GF2_Vector *out, const BPU_T_GF2_Vector *in, const struct _BPU_T_Code_Ctx *ctx) {
	int i, j, bit_in_msg = 0;
	int ret = 0;

	BPU_T_GF2_Poly *code_msg_temp;
	BPU_T_GF2_Poly code_msg;

	ret += BPU_gf2PolyMalloc(&code_msg, ctx->code_spec->qcldpc->m);
	code_msg_temp = (BPU_T_GF2_Poly*)calloc(ctx->code_spec->qcldpc->n0, sizeof(BPU_T_GF2_Poly));

	for (i = 0; i < ctx->code_spec->qcldpc->n0; i++) 
		ret += BPU_gf2PolyMalloc(&code_msg_temp[i], ctx->code_spec->qcldpc->m);
	if (ret != 0) {
		BPU_printError("Could not allocate memory");
		return -1;
	}
	
	BPU_printDebug("calculating plaintext x G_masked");
	//for every row of blocks in G_masked
	for (i = 0; i < ctx->code_spec->qcldpc->n0 - 1; i++) {
		//for every row in a block
		for (j = 0; j < ctx->code_spec->qcldpc->m; j++) {
			if (BPU_gf2PolyGetBit(in, bit_in_msg)) {
				//xor row from the main diagonal block in G_masked
				BPU_gf2VecXor(&code_msg_temp[i], &ctx->code_spec->qcldpc->G_masked.matrices[i]);
				//xor row from the blocks not on the main diagonal in G_masked
				BPU_gf2VecXor(&code_msg_temp[ctx->code_spec->qcldpc->n0 - 1], &ctx->code_spec->qcldpc->G_masked.matrices[i + (ctx->code_spec->qcldpc->n0 - 1)]);
			}
			//cyclic shift by 1
			BPU_gf2PolyMulX(&ctx->code_spec->qcldpc->G_masked.matrices[i]);
			BPU_gf2PolyMulX(&ctx->code_spec->qcldpc->G_masked.matrices[i + (ctx->code_spec->qcldpc->n0 - 1)]);

			bit_in_msg++;
		}
	}

	//concatenate partial code message blocks
	for (i = ctx->code_spec->qcldpc->n0 - 1; i >= 0; i--)	{
		BPU_gf2PolyAdd(&code_msg, &code_msg_temp[i], 0);
		BPU_gf2PolyFree(&code_msg_temp[i], 0);
		if (i > 0)
			BPU_gf2PolyShiftLeft(&code_msg, ctx->code_spec->qcldpc->m);
	}

	BPU_gf2PolyAdd(out, &code_msg, 0);
	
#if defined(DEBUG_L)
	BPU_printDebug("code message: ");
	//BPU_printGf2Poly(out);
#endif

	BPU_printDebug("free and exit");
	free(code_msg_temp);
	BPU_gf2PolyFree(&code_msg, 0);

	return 0;
}

int BPU_mecsQcldpcEncode2(BPU_T_GF2_Vector *out, const BPU_T_GF2_Vector *in, const struct _BPU_T_Code_Ctx *ctx) {
	int32_t i, bit, row, block, bit_in_msg = 0;
	int32_t ret = 0;

	BPU_T_GF2_Poly *code_msg_temp;
	BPU_T_GF2_Poly code_msg;

	ret += BPU_gf2PolyMalloc(&code_msg, ctx->code_spec->qcldpc->m);
	code_msg_temp = (BPU_T_GF2_Poly*)calloc(ctx->code_spec->qcldpc->n0, sizeof(BPU_T_GF2_Poly));

	for (i = 0; i < ctx->code_spec->qcldpc->n0; i++)
		ret += BPU_gf2PolyMalloc(&code_msg_temp[i], ctx->code_spec->qcldpc->m);
	if (ret != 0) {
		BPU_printError("Could not allocate memory");
		return -1;
	}
	
	BPU_printDebug("calculating plaintext x G_masked");
	//for every row of blocks in G_masked
	for (row = 0; row < ctx->code_spec->qcldpc->n0 - 1; row++) {
		//check every bit for in plaintext for 1
		for (bit = 0; bit < ctx->code_spec->qcldpc->m; bit++) {
			//if bit==1, xor corresponding row in G_masked and code_message
			if (BPU_gf2PolyGetBit(in, bit_in_msg))
			{
				for (block = 0; block < ctx->code_spec->qcldpc->n0; block++)
					BPU_gf2VecXor(&code_msg_temp[block], &ctx->code_spec->qcldpc->G_masked.matrices[row * ctx->code_spec->qcldpc->n0 + block]);
			}
			//shift row in G_masked
			for (block = 0; block < ctx->code_spec->qcldpc->n0; block++)
				BPU_gf2PolyMulX(&ctx->code_spec->qcldpc->G_masked.matrices[row * ctx->code_spec->qcldpc->n0 + block]);
			bit_in_msg++;
		}
	}

	//concatenate partial code message blocks
	for (i = ctx->code_spec->qcldpc->n0 - 1; i >= 0; i--)	{
		BPU_gf2PolyAdd(&code_msg, &code_msg_temp[i], 0);
		BPU_gf2PolyFree(&code_msg_temp[i], 0);
		if (i > 0)
			BPU_gf2PolyShiftLeft(&code_msg, ctx->code_spec->qcldpc->m);
	}
	
	BPU_gf2PolyAdd(out, &code_msg, 0);
	
#if defined(DEBUG_L)
	BPU_printDebug("code message: ");
	//BPU_printGf2Poly(out);
#endif
	
	BPU_printDebug("free and exit");
	free(code_msg_temp);
	BPU_gf2PolyFree(&code_msg, 0);

	return 0;
}
#endif

#ifdef BPU_CONF_DECRYPTION
int BPU_mecsQcldpcDecrypt1(BPU_T_GF2_Vector *out, const BPU_T_GF2_Vector *in, const struct _BPU_T_Code_Ctx *ctx) {

	int32_t ret = 0, delta = BPU_QCLDPC_PARAM_DELTA;
	int32_t i, bit;
	int32_t bit_in_msg = 0;

	BPU_T_GF2_Poly *ct_Q_temp, *ct_Q_S_temp;
	BPU_T_GF2_Poly ct_Q, ct_Q_S;
	BPU_T_GF2_Sparse_Poly row;

	ct_Q_temp = (BPU_T_GF2_Poly*)calloc(ctx->code_spec->qcldpc->n0, sizeof(BPU_T_GF2_Poly));
	ct_Q_S_temp = (BPU_T_GF2_Poly*)calloc(ctx->code_spec->qcldpc->n0 - 1, sizeof(BPU_T_GF2_Poly));

	ret += BPU_gf2PolyMalloc(&ct_Q_S, ctx->code_spec->qcldpc->m);
	ret += BPU_gf2PolyMalloc(&ct_Q, ctx->code_spec->qcldpc->m);
	for(i=0;i<ctx->code_spec->qcldpc->n0;i++)
		ret += BPU_gf2PolyMalloc(&ct_Q_temp[i], ctx->code_spec->qcldpc->m);
	if (ret != 0) {
		BPU_printError("Could not allocate memory");
		return -1;
	}
	
	BPU_printDebug("calculating ciphertext x Q");
	//for every row of blocks in Q
	for (i = 0; i < ctx->code_spec->qcldpc->n0; i++) {		
		//check every bit in ciphertext for 1
		for (bit = 0; bit < ctx->code_spec->qcldpc->m; bit++) {
			//if bit==1, xor corresponding row in Q
			if (BPU_gf2PolyGetBit(in, bit_in_msg)) {				
				BPU_gf2SparseQcMatrixGetRow(&row, &ctx->code_spec->qcldpc->Q_sparse, bit_in_msg);
				BPU_gf2SparsePolyAdd(&ct_Q_temp[i], &row);
				BPU_gf2SparsePolyFree(&row, 0);
			}
			bit_in_msg++;
		}
	}

	//concatenate partial blocks
	for (i = ctx->code_spec->qcldpc->n0 - 1; i >= 0; i--)
	{
		BPU_gf2PolyAdd(&ct_Q, &ct_Q_temp[i], 0);
		BPU_gf2PolyFree(&ct_Q_temp[i], 0);
		if (i > 0)
			BPU_gf2PolyShiftLeft(&ct_Q, ctx->code_spec->qcldpc->m);
	}
	
#if defined(DEBUG_L)
	BPU_printDebug("ciphertext x Q: ");
	//BPU_printGf2Poly(&ct_Q);
#endif

	// null error vector
	BPU_gf2VecNull(ctx->e);

	// try to decode with faster algorithm
	if (!BPU_mecsQcldpcDecode2(ctx->e, &ct_Q, ctx)) {
		// free decoded
		BPU_gf2VecNull(ctx->e);
		while (1) {
			// if not decoded, try algorithm with lower DFR
			if (!BPU_mecsQcldpcDecode1(ctx->e, &ct_Q, delta, ctx)) {
				// free decoded
				BPU_gf2VecNull(ctx->e);
				// if not decoded decrease threshold tolerance param
				delta--;
				if (delta < 0) {
					BPU_printError("Decoding failed; delta is zero");
					return -1;						
				}
			}
			else
				break;
		}
	}
		
	BPU_gf2VecXor(&ct_Q, ctx->e);
	BPU_gf2VecCrop(out, &ct_Q, 0, ctx->msg_len);
	
	bit_in_msg = 0;
	BPU_printDebug("calculating decoded_message x S");
	//for every row of blocks in S
	for (i = 0; i < ctx->code_spec->qcldpc->n0 - 1; i++) {
		BPU_gf2PolyMalloc(&ct_Q_S_temp[i], ctx->code_spec->qcldpc->m);
		//check every bit in decoded_message for 1
		for (bit = 0; bit < ctx->code_spec->qcldpc->m; bit++) {
			//if bit==1, xor corresponding row in S
			if (BPU_gf2PolyGetBit(out, bit_in_msg)) {				
				BPU_gf2VecXor(&ct_Q_S_temp[i], &ctx->code_spec->qcldpc->S.matrices[i]);
			}
			//cyclic shift by 1
			BPU_gf2PolyMulX(&ctx->code_spec->qcldpc->S.matrices[i]);

			bit_in_msg++;
		}
	}

	//concatenate partial blocks
	for (i = ctx->code_spec->qcldpc->n0 - 2; i >= 0; i--)
	{
		BPU_gf2PolyAdd(&ct_Q_S, &ct_Q_S_temp[i], 0);
		BPU_gf2PolyFree(&ct_Q_S_temp[i], 0);
		if (i > 0)
			BPU_gf2PolyShiftLeft(&ct_Q_S, ctx->code_spec->qcldpc->m);
	}

	BPU_gf2VecNull(out);
	BPU_gf2PolyAdd(out, &ct_Q_S, 0);

#if defined(DEBUG_L)
	BPU_printDebug("decoded_message x S: ");
	//BPU_printGf2Poly(out);
#endif
	
	BPU_printDebug("free and exit");
	free(ct_Q_S_temp);
	free(ct_Q_temp);
	BPU_gf2PolyFree(&ct_Q, 0);
	BPU_gf2PolyFree(&ct_Q_S, 0);

	return 0;
}

int BPU_mecsQcldpcDecrypt2(BPU_T_GF2_Vector *out, const BPU_T_GF2_Vector *in, const struct _BPU_T_Code_Ctx *ctx) {
	// printf("BPU_mecsQcldpcDecrypt2\n");

	int32_t ret = 0, delta = BPU_QCLDPC_PARAM_DELTA;
	int32_t i;
	int32_t bit_in_msg = 0;
	int row, bit,block ;

	BPU_T_GF2_Poly *ct_Q_temp, *ct_Q_S_temp;
	BPU_T_GF2_Poly ct_Q, ct_Q_S;
	BPU_T_GF2_Sparse_Poly row_poly;

	ct_Q_temp = (BPU_T_GF2_Poly*)calloc(ctx->code_spec->qcldpc->n0, sizeof(BPU_T_GF2_Poly));
	ct_Q_S_temp = (BPU_T_GF2_Poly*)calloc(ctx->code_spec->qcldpc->n0 - 1, sizeof(BPU_T_GF2_Poly));

	ret += BPU_gf2PolyMalloc(&ct_Q_S, ctx->code_spec->qcldpc->m);
	ret += BPU_gf2PolyMalloc(&ct_Q, ctx->code_spec->qcldpc->m);
	
	for (i = 0; i < ctx->code_spec->qcldpc->n0; i++)
		ret += BPU_gf2PolyMalloc(&ct_Q_temp[i], ctx->code_spec->qcldpc->m);
	for (i = 0; i < ctx->code_spec->qcldpc->n0 -1 ; i++)
		ret += BPU_gf2PolyMalloc(&ct_Q_S_temp[i], ctx->code_spec->qcldpc->m);
	
	if (ret != 0) {
		BPU_printError("Could not allocate memory");
		return -1;
	}
	
	BPU_printDebug("calculating ciphertext x Q");
	for (row = 0; row < ctx->code_spec->qcldpc->n0; row++) {
		for (bit = 0; bit < ctx->code_spec->qcldpc->m; bit++) {
			if (BPU_gf2PolyGetBit(in, bit_in_msg))
				for (block = 0; block < ctx->code_spec->qcldpc->n0; block++) {
					BPU_gf2SparsePolyGetShift(&row_poly, &ctx->code_spec->qcldpc->Q_sparse.matrices[block + row * ctx->code_spec->qcldpc->n0], bit, ctx->code_spec->qcldpc->m);
					BPU_gf2SparsePolyAdd(&ct_Q_temp[block], &row_poly);
					BPU_gf2SparsePolyFree(&row_poly, 0);
				}

			bit_in_msg++;
		}
	}

	//concatenate partial blocks
	for (i = ctx->code_spec->qcldpc->n0 - 1; i >= 0; i--) {
		BPU_gf2PolyAdd(&ct_Q, &ct_Q_temp[i], 0);
		BPU_gf2PolyFree(&ct_Q_temp[i], 0);
		if (i > 0)
			BPU_gf2PolyShiftLeft(&ct_Q, ctx->code_spec->qcldpc->m);
	}
		
#if defined(DEBUG_L)
	BPU_printDebug("ciphertext x Q: ");
	//BPU_printGf2Poly(&ct_Q);
#endif

	// null error vector
	BPU_gf2VecNull(ctx->e);

	// try to decode with faster algorithm
	if (!BPU_mecsQcldpcDecode2(ctx->e, &ct_Q, ctx)) {
		// free decoded
		BPU_gf2VecNull(ctx->e);
		while (1) {
			// if not decoded, try algorithm with lower DFR
			if (!BPU_mecsQcldpcDecode1(ctx->e, &ct_Q, delta, ctx)) {
				// free decoded
				BPU_gf2VecNull(ctx->e);
				// if not decoded decrease threshold tolerance param
				delta--;
				if (delta < 0) {
					// BPU_printError("Decoding failed; delta is zero");
					break;
					//return -1;					
				}
			}
			else
				break;
		}
	}
	
	BPU_gf2VecXor(&ct_Q, ctx->e);
	BPU_gf2VecCrop(out, &ct_Q, 0, ctx->msg_len);
	
	bit_in_msg = 0;
	BPU_printDebug("calculating decoded_message x S");
	
	for (row = 0; row < ctx->code_spec->qcldpc->n0 - 1; row++) {
		for (bit = 0; bit < ctx->code_spec->qcldpc->m; bit++) {
			if (BPU_gf2PolyGetBit(out, bit_in_msg))
				for (block = 0; block < ctx->code_spec->qcldpc->n0 - 1; block++)
					BPU_gf2VecXor(&ct_Q_S_temp[block], &ctx->code_spec->qcldpc->S.matrices[block + row * (ctx->code_spec->qcldpc->n0 - 1)]);
	
			for (block = 0; block < ctx->code_spec->qcldpc->n0 - 1; block++)
				BPU_gf2PolyMulX(&ctx->code_spec->qcldpc->S.matrices[block + row * (ctx->code_spec->qcldpc->n0 - 1)]);
	
			bit_in_msg++;
		}
	}

	//concatenate partial blocks
	for (i = ctx->code_spec->qcldpc->n0 - 2; i >= 0; i--) {
		BPU_gf2PolyAdd(&ct_Q_S, &ct_Q_S_temp[i], 0);
		BPU_gf2PolyFree(&ct_Q_S_temp[i], 0);
		if (i > 0)
			BPU_gf2PolyShiftLeft(&ct_Q_S, ctx->code_spec->qcldpc->m);
	}

	BPU_gf2VecNull(out);
	BPU_gf2PolyAdd(out, &ct_Q_S, 0);

#if defined(DEBUG_L)
	BPU_printDebug("decoded_message x S: ");
	//BPU_printGf2Poly(out);
#endif
	
	BPU_printDebug("free and exit");
	free(ct_Q_S_temp);
	free(ct_Q_temp);
	BPU_gf2PolyFree(&ct_Q, 0);
	BPU_gf2PolyFree(&ct_Q_S, 0);

	return 0;
}

int BPU_mecsQcldpcDecode1(BPU_T_GF2_Vector *error_vec, const BPU_T_GF2_Vector *cipher_text, int delta, const struct _BPU_T_Code_Ctx *ctx) {
	// printf("BPU_mecsQcldpcDecode1\n");

	BPU_T_GF2_Poly syndrom;
	BPU_T_GF2_Sparse_Poly row;
	int32_t iter = -1, max, upc, isSyndromZero = 0;
	int32_t flipped_bits = 0, flipped_bits_iter = 0;
	uint8_t *upc_counts; //TODO optimize memory to use uint16_t
	int32_t bit;

	upc_counts = (uint8_t*)calloc(cipher_text->len, sizeof(uint8_t));

	// calc the syndrom
	BPU_mecsQcldpcCalcSyndrom(&syndrom, cipher_text, ctx);
	// check syndrom
	if (!BPU_gf2PolyIsZero(&syndrom)) {
		// for max iterations
		for (iter = 0; iter < BPU_QCLDPC_PARAM_MAX_ITER; iter++) {
			max = 0;
			// for every bit of cipher text
			for (bit = 0; bit < error_vec->len; bit++) {
				// calc #UPC
				BPU_gf2SparseQcMatrixGetRow(&row, &ctx->code_spec->qcldpc->H, bit);
				upc = BPU_gf2SparsePolyAndHW(&syndrom, &row);
				upc_counts[bit] = upc;
				if (upc > max)
					max = upc;
				BPU_gf2SparsePolyFree(&row, 0);
			}

			if (max == 0) {
				isSyndromZero = 0;
				break;
			}

			flipped_bits_iter = 0;
			// check which bits to flip
			for (bit = 0; bit < error_vec->len; bit++) {
				if (upc_counts[bit] > 0 && upc_counts[bit] >= (max - delta)) {
					flipped_bits++; flipped_bits_iter++;
					// flip bit
					BPU_gf2VecSetBit(error_vec, bit, !BPU_gf2VecGetBit(error_vec, bit));
					// update syndrom
					BPU_gf2SparseQcMatrixGetRow(&row, &ctx->code_spec->qcldpc->H, bit);
					BPU_gf2SparsePolyAdd(&syndrom, &row);
					BPU_gf2SparsePolyFree(&row, 0);
					// check the syndrom
					if (BPU_gf2PolyIsZero(&syndrom)) {
						isSyndromZero = 1;
						break;
					}
				}
			}

			if (isSyndromZero)
				break;
		}
	}
	else {
		isSyndromZero = 1;
	}
	
	//free
	BPU_gf2PolyFree(&syndrom, 0);
	free(upc_counts);

	return isSyndromZero;
}

int BPU_mecsQcldpcDecode2(BPU_T_GF2_Vector *error_vec, const BPU_T_GF2_Vector *cipher_text, const struct _BPU_T_Code_Ctx *ctx) {
	// printf("BPU_mecsQcldpcDecode2\n");
	
	BPU_T_GF2_Poly syndrom;
	BPU_T_GF2_Sparse_Poly row;
	int32_t iter = -1, upc, isSyndromZero = 0;
	int32_t flipped_bits = 0, flipped_bits_iter = 0;
	const uint16_t B_store[BPU_QCLDPC_MAX_B_VALUES] = { 28,26,24,22,20 };

	int32_t bit;
	// calc the syndrom
	BPU_mecsQcldpcCalcSyndrom(&syndrom, cipher_text, ctx);

	// check syndrom
	if (!BPU_gf2PolyIsZero(&syndrom)) {
		// for max iterations
		for (iter = 0; iter < BPU_QCLDPC_PARAM_MAX_ITER; iter++) {
			flipped_bits_iter = 0;
			// for every bit of cipher text
			for (bit = 0; bit < error_vec->len; bit++) {
				// calc #UPC
				BPU_gf2SparseQcMatrixGetRow(&row, &ctx->code_spec->qcldpc->H, bit);
				upc = BPU_gf2SparsePolyAndHW(&syndrom, &row);

				// check which bits to flip
				if (upc > 0 && upc >= B_store[iter < BPU_QCLDPC_MAX_B_VALUES ? iter : (BPU_QCLDPC_MAX_B_VALUES - 1)] - BPU_QCLDPC_PARAM_DELTA_B) {
					flipped_bits++; flipped_bits_iter++;
					// flip bit
					BPU_gf2VecSetBit(error_vec, bit, !BPU_gf2VecGetBit(error_vec, bit));
					// update syndrom
					BPU_gf2SparsePolyAdd(&syndrom, &row);
					// check the syndrom
					if (BPU_gf2PolyIsZero(&syndrom)) {
						isSyndromZero = 1;
						BPU_gf2SparsePolyFree(&row, 0);
						break;
					}
				}
				BPU_gf2SparsePolyFree(&row, 0);
			}

			if (flipped_bits_iter < 1) {
				isSyndromZero = 0;
				break;
			}

			if (isSyndromZero)
				break;
		}
	}
	else {
		isSyndromZero = 1;
	}
	
	//free
	BPU_gf2PolyFree(&syndrom, 0);

	return isSyndromZero;
}

void BPU_mecsQcldpcCalcSyndrom(BPU_T_GF2_Vector *syndrom, const BPU_T_GF2_Vector *cipher_text, const struct _BPU_T_Code_Ctx *ctx) {
	BPU_T_GF2_Sparse_Poly row;
	int32_t i;

	BPU_gf2PolyMalloc(syndrom, ctx->code_spec->qcldpc->H.n);
	for (i = 0; i < cipher_text->len; i++) {
		if (BPU_gf2VecGetBit(cipher_text, i) == 1ul) {
			BPU_gf2SparseQcMatrixGetRow(&row, &ctx->code_spec->qcldpc->H, i);
			BPU_gf2SparsePolyAdd(syndrom, &row);
			BPU_gf2SparsePolyFree(&row, 0);
		}
	}
}
#endif

#ifdef BPU_CONF_KEY_GEN
int BPU_mecsQcldpcGenKeys(BPU_T_Code_Ctx *ctx) {
	
	BPU_T_GF2_Poly mod;
	BPU_T_GF2_Poly H_last_inv, test_inv;
	BPU_T_GF2_Poly poly_temp, poly_inv_temp;

	BPU_T_GF2_Poly *G_temp, *H_temp;
	BPU_T_GF2_QC_Matrix G_temp_mat, H_temp_mat;
	BPU_T_GF2_Sparse_Qc_Matrix H_temp_sparse;

	int32_t i, ret = 0, err;
	int32_t *wi = (int*)calloc(ctx->code_spec->qcldpc->n0, sizeof(int32_t));

	G_temp = (BPU_T_GF2_Poly*)calloc(ctx->code_spec->qcldpc->n0 - 1, sizeof(BPU_T_GF2_Poly));
	H_temp = (BPU_T_GF2_Poly*)calloc(ctx->code_spec->qcldpc->n0, sizeof(BPU_T_GF2_Poly));
	
	BPU_gf2PolyMalloc(&H_last_inv, ctx->code_spec->qcldpc->m);

	// init modulus like 1000000...0001
	BPU_gf2PolyMalloc(&mod, ctx->code_spec->qcldpc->m + 1);
	BPU_gf2VecSetBit(&mod, ctx->code_spec->qcldpc->m, 1ul);
	BPU_gf2VecSetBit(&mod, 0, 1ul);


#if defined(DEBUG_L)
	BPU_printDebug("modulus: ");
	//BPU_printGf2Poly(&mod);
	BPU_printDebug("generating H vectors");
#endif

	// alloc parity-check matrix
	err = 0;
	for (i = 0; i < ctx->code_spec->qcldpc->n0; i++) {
		// calc weight of polynomials (last poly must be odd)
		if ((ctx->code_spec->qcldpc->w / ctx->code_spec->qcldpc->n0) % 2 == 1)
			wi[i] = ctx->code_spec->qcldpc->w / ctx->code_spec->qcldpc->n0 + (int)(i < (ctx->code_spec->qcldpc->w%ctx->code_spec->qcldpc->n0));
		else
			wi[i] = ctx->code_spec->qcldpc->w / ctx->code_spec->qcldpc->n0 + (int)(i < (ctx->code_spec->qcldpc->w%ctx->code_spec->qcldpc->n0)) + (int)(i == 0) - (int)(i == ctx->code_spec->qcldpc->n0 - 1);

		// generate random polynomials of given weight	
		err += BPU_gf2PolyInitRand(&H_temp[i], ctx->code_spec->qcldpc->m, wi[i], 0);

#if defined(DEBUG_L)
		BPU_printDebug("H[%i]: ", i);
		//BPU_printGf2Poly(&H_temp[i]);
#endif
	}

	if (!err)
		BPU_printDebug("generation successful");
	else
		BPU_printError("generation failed");


	BPU_printDebug("finding inversion to H[%i]", ctx->code_spec->qcldpc->n0 - 1);
	// check if H[n0-1] has inversion
	while (!ret) {

		// find inversion using XGCD
		BPU_gf2PolyCopy(&poly_temp, &H_temp[ctx->code_spec->qcldpc->n0 - 1]);
		BPU_gf2PolySetDeg(&poly_temp, -1);
		ret = BPU_gf2PolyInv(&poly_inv_temp, &poly_temp, &mod);

		// if inversion exists, test it (poly x inversion modulo = 1)
		if (ret) {

			BPU_printDebug("testing inversion");
			BPU_gf2PolyMulMod(&poly_inv_temp, &poly_temp, &test_inv, &mod, 1);

			if (test_inv.len != 1 || test_inv.elements[0] != 1ul) {
				ret = 0;
				BPU_printWarning("inversion failed");
			}
			else {
				BPU_printDebug("inversion OK");
				BPU_gf2PolyAdd(&H_last_inv, &poly_inv_temp, 0);
				BPU_gf2PolyFree(&poly_temp, 0);
				BPU_gf2PolyFree(&poly_inv_temp, 0);
			}

			BPU_gf2PolyFree(&test_inv, 0);
		}

		// inversion not found, regenerate last poly and try to find inversion again
		if (!ret) {
			BPU_printDebug("inversion not found");
			BPU_printDebug("generating new H[%i]", ctx->code_spec->qcldpc->n0 - 1);
			BPU_gf2PolyFree(&H_temp[ctx->code_spec->qcldpc->n0 - 1], 0);
			ret += BPU_gf2PolyInitRand(&H_temp[ctx->code_spec->qcldpc->n0 - 1], ctx->code_spec->qcldpc->m, wi[ctx->code_spec->qcldpc->n0 - 1], 0);

#if defined(DEBUG_L)
			BPU_printDebug("H[%i]: ", ctx->code_spec->qcldpc->n0 - 1);
			//BPU_printGf2Poly(&H_temp[ctx->code_spec->qcldpc->n0-1]);
#endif
		}
	}
#if defined(DEBUG_L)
	BPU_printDebug("inversion to H[%i] found: ", ctx->code_spec->qcldpc->n0 - 1);
	//BPU_printGf2Poly(&H_last_inv);
	BPU_printDebug("creating H matrix");
#endif

	// create H temp matrix
	BPU_gf2QcMatrixMalloc(&H_temp_mat, ctx->code_spec->qcldpc->n0, ctx->code_spec->qcldpc->m, 1, 0);
	for (i = 0; i < ctx->code_spec->qcldpc->n0; i++) {
		BPU_gf2PolyCopy(&H_temp_mat.matrices[i], &H_temp[i]);
	}
	BPU_gf2QcMatrixToSparse(&H_temp_sparse, &H_temp_mat, wi);	

	// create G temp matrix
	for (i = 0; i < ctx->code_spec->qcldpc->n0 - 1; i++) {
		BPU_printDebug("multiplicating vectors H[%i]^-1 x H[%i]", ctx->code_spec->qcldpc->n0 - 1, i);
		BPU_gf2PolyMulMod(&H_last_inv, &H_temp[i], &G_temp[i], &mod, 0);
		BPU_gf2PolyFree(&H_temp[i], 0);
#if defined(DEBUG_L)
		BPU_printDebug("G[%i]: ", i);
		//BPU_printGf2Poly(&G_temp[i]);
#endif
	}
	BPU_gf2PolyFree(&H_temp[ctx->code_spec->qcldpc->n0 - 1], 0);
	BPU_gf2PolyFree(&H_last_inv, 0);

	BPU_printDebug("creating temp G for GH^T test");
	BPU_gf2QcMatrixMalloc(&G_temp_mat, ctx->code_spec->qcldpc->n0 - 1, ctx->code_spec->qcldpc->m, 0, 1);
	for (i = 0; i < ctx->code_spec->qcldpc->n0 - 1; i++) {
		BPU_gf2PolyCopy(&G_temp_mat.matrices[i], &G_temp[i]);
	}
	
	// fprintf(stderr, "BPU_mecsQcldpcGenKeys\n");
	// fprintf(stderr, " - H.n: %i\n", H_temp_sparse.n);
	// fprintf(stderr, " - H.k: %i\n", H_temp_sparse.k);
	// fprintf(stderr, " - G.n: %i\n", G_temp_mat.n);
	// fprintf(stderr, " - G.k: %i\n\n", G_temp_mat.k);

	BPU_printDebug("testing GH^T");
	ret = 0;
	// test if G x H^T = 0
	if (BPU_mecsQcldpcTestGHmatrices(&G_temp_mat, &H_temp_sparse) != 0) {
		BPU_printError("generator x parity check matrix ERROR");
		ret = -1;
	}
	else {
		BPU_printDebug("GH^t = 0");
	}
	
	BPU_printDebug("transposing G polynomials");
	if (ret == 0) {
		for (i = 0; i < ctx->code_spec->qcldpc->n0 - 1; i++) {
			BPU_gf2PolyFree(&G_temp[i], 0);
			BPU_gf2PolyTransp(&G_temp[i], &G_temp_mat.matrices[i]);

#if defined(DEBUG_L)
			BPU_printDebug("G[%i]: ", i);
			//BPU_printGf2Poly(&G_temp[i]);
#endif
		}
	}
	
	BPU_printDebug("transposing H matrix");
	if (ret == 0)
		BPU_gf2SparseQcMatrixTransp(&ctx->code_spec->qcldpc->H, &H_temp_sparse);

	free(H_temp);
	free(wi);
	BPU_gf2PolyFree(&mod, 0);
	BPU_gf2QcMatrixFree(&G_temp_mat, 0);
	BPU_gf2QcMatrixFree(&H_temp_mat, 0);
	BPU_gf2SparseQcMatrixFree(&H_temp_sparse, 0);

	BPU_printDebug("generating masked G");
	if (ret == 0) {
		if (ctx->type == BPU_EN_CODE_QCLDPC1) { 
			ret += BPU_QcldpcMaskKeys1(ctx, G_temp);
		}
		if (ctx->type == BPU_EN_CODE_QCLDPC2) {
			ret += BPU_QcldpcMaskKeys2(ctx, G_temp);  
		}
	}
	BPU_printDebug("free and exit");

	free(G_temp);
	return ret;
}


int BPU_QcldpcMaskKeys1(BPU_T_Code_Ctx *ctx, BPU_T_GF2_Poly *G_temp) {
	BPU_T_GF2_Poly mod;
	BPU_T_GF2_Poly test_inv;
	BPU_T_GF2_Poly poly_temp, poly_inv_temp;
	BPU_T_GF2_Poly *S_temp, *S_inv_temp, *Q_temp, *Q_inv_temp;
	BPU_T_GF2_Poly *SG_temp, *G_masked;
	BPU_T_GF2_QC_Matrix Q_temp_mat;

	S_temp = (BPU_T_GF2_Poly*)calloc(ctx->code_spec->qcldpc->n0 - 1, sizeof(BPU_T_GF2_Poly));
	S_inv_temp = (BPU_T_GF2_Poly*)calloc(ctx->code_spec->qcldpc->n0 - 1, sizeof(BPU_T_GF2_Poly));
	Q_temp = (BPU_T_GF2_Poly*)calloc(ctx->code_spec->qcldpc->n0, sizeof(BPU_T_GF2_Poly));
	Q_inv_temp = (BPU_T_GF2_Poly*)calloc(ctx->code_spec->qcldpc->n0, sizeof(BPU_T_GF2_Poly));
	SG_temp = (BPU_T_GF2_Poly*)calloc((ctx->code_spec->qcldpc->n0 - 1) * 2, sizeof(BPU_T_GF2_Poly));
	G_masked = (BPU_T_GF2_Poly*)calloc((ctx->code_spec->qcldpc->n0 - 1) * 2, sizeof(BPU_T_GF2_Poly));

	int32_t i, j, ret, err;
	//Hamming weight of Matrix Q blocks
	int32_t wq[] = { 11,11,11 };

	// init modulus like 1000000...0001
	BPU_gf2PolyMalloc(&mod, ctx->code_spec->qcldpc->m + 1);
	BPU_gf2VecSetBit(&mod, ctx->code_spec->qcldpc->m, 1ul);
	BPU_gf2VecSetBit(&mod, 0, 1ul);

#if defined(DEBUG_L)
	BPU_printDebug("modulus: ");
	//BPU_printGf2Poly(&mod);	
#endif

	err = 0;
	//Generate matrix S
	//Generating only cyclic blocks along the main diagonal, elsewhere are null matrix blocks
	BPU_printDebug("generating matrix S polynomials");
	for (i = 0; i < ctx->code_spec->qcldpc->n0 - 1; i++)
	{
		err += BPU_gf2PolyInitRandProb(&S_temp[i], ctx->code_spec->qcldpc->m, 0.5, 0, 0);
		if (err)
		{
			BPU_printError("S generation failed");
			return -1;
		}

#if defined(DEBUG_L)
		BPU_printDebug("S[%i]: ", i);
		//BPU_printGf2Poly(&S_temp[i]);
#endif
	}


	BPU_printDebug("inverting matrix S");
	for (i = 0; i < ctx->code_spec->qcldpc->n0 - 1; i++)	{
		ret = 0;
		BPU_gf2PolyMalloc((&S_inv_temp[i]), ctx->code_spec->qcldpc->m);
		while (!ret)	{
			BPU_gf2PolyCopy(&poly_temp, &S_temp[i]);
			BPU_gf2PolySetDeg(&poly_temp, -1);
			ret = BPU_gf2PolyInv(&poly_inv_temp, &poly_temp, &mod);
			if (ret)	{
				BPU_gf2PolyMulMod(&poly_inv_temp, &poly_temp, &test_inv, &mod, 1);
				if (test_inv.len != 1 || test_inv.elements[0] != 1ul) {
					ret = 0;
					BPU_printWarning("inversion failed");
				}
				else {
					BPU_printDebug("inversion S[%i] OK",i);
					BPU_gf2PolyAdd(&S_inv_temp[i], &poly_inv_temp, 0);
					BPU_gf2PolyFree(&poly_inv_temp, 0);
					BPU_gf2PolyFree(&poly_temp, 0);
				}
				BPU_gf2PolyFree(&test_inv, 0);
			}

			if (!ret)	{
				BPU_printDebug("inversion S[%i] not found",i);
				BPU_printDebug("generating new S[%i]", i);
				BPU_gf2PolyFree(&poly_inv_temp, 0);
				BPU_gf2PolyFree(&poly_temp, 0);
				BPU_gf2PolyFree(&S_temp[i], 0);
				ret += BPU_gf2PolyInitRandProb(&S_temp[i], ctx->code_spec->qcldpc->m, 0.5, 0, 0);

#if defined(DEBUG_L)
				BPU_printDebug("S[%i]: ", i);
				//BPU_printGf2Poly(&S_temp[i]);
#endif
			}
		}
	}

	BPU_printDebug("copying matrix S into code context");
	BPU_gf2QcMatrixMalloc(&ctx->code_spec->qcldpc->S, ctx->code_spec->qcldpc->n0 - 1, ctx->code_spec->qcldpc->m, 0, 0);
	for (i = 0; i < ctx->code_spec->qcldpc->n0 - 1; i++) {
		BPU_gf2PolyCopy(&ctx->code_spec->qcldpc->S.matrices[i], &S_temp[i]);
		BPU_gf2PolyFree(&S_temp[i], 0);
	}

	//Generating matrix Q
	//Generating only cyclic blocks along the main diagonal, elsewhere are zero matrix blocks
	err = 0;
	BPU_printDebug("generating matrix Q polynomials");
	for (i = 0; i < ctx->code_spec->qcldpc->n0; i++)	{
		err += BPU_gf2PolyInitRand(&Q_temp[i], ctx->code_spec->qcldpc->m, wq[i], 0);
		if (err)	{
			BPU_printError("Q generation failed");
			return -1;
		}

#if defined(DEBUG_L)
		BPU_printDebug("Q[%d]: ", i);
		//BPU_printGf2Poly(&Q_temp[i]);
#endif
	}


	BPU_printDebug("inverting matrix Q");
	for (i = 0; i < ctx->code_spec->qcldpc->n0; i++)	{
		ret = 0;
		BPU_gf2PolyMalloc((&Q_inv_temp[i]), ctx->code_spec->qcldpc->m);
		while (!ret)	{
			BPU_gf2PolyCopy(&poly_temp, &Q_temp[i]);
			BPU_gf2PolySetDeg(&poly_temp, -1);
			ret = BPU_gf2PolyInv(&poly_inv_temp, &poly_temp, &mod);
			if (ret)	{
				BPU_gf2PolyMulMod(&poly_inv_temp, &poly_temp, &test_inv, &mod, 1);
				if (test_inv.len != 1 || test_inv.elements[0] != 1ul) {
					ret = 0;
					BPU_printWarning("inversion failed");
				}
				else {
					BPU_printDebug("inversion Q[%i] OK",i);
					BPU_gf2PolyAdd(&Q_inv_temp[i], &poly_inv_temp, 0);
					BPU_gf2PolyFree(&poly_inv_temp, 0);
					BPU_gf2PolyFree(&poly_temp, 0);
				}
				BPU_gf2PolyFree(&test_inv, 0);
			}

			if (!ret)	{
				BPU_printDebug("inversion Q[%i] not found",i);
				BPU_printDebug("generating new Q[%i]", i);
				BPU_gf2PolyFree(&poly_inv_temp, 0);
				BPU_gf2PolyFree(&poly_temp, 0);
				BPU_gf2PolyFree(&Q_temp[i], 0);
				ret += BPU_gf2PolyInitRand(&Q_temp[i], ctx->code_spec->qcldpc->m, wq[i], 0);

#if defined(DEBUG_L)
				BPU_printDebug("Q[%i]: ", i);
				//BPU_printGf2Poly(&Q_temp[i]);
#endif
			}
		}
	}
	
	BPU_printDebug("copying sparse matrix Q into code context");
	BPU_gf2QcMatrixMalloc(&Q_temp_mat, ctx->code_spec->qcldpc->n0, ctx->code_spec->qcldpc->m, 0, 0);
	for (i = 0; i < ctx->code_spec->qcldpc->n0; i++) {
		BPU_gf2PolyCopy(&Q_temp_mat.matrices[i], &Q_temp[i]);
		BPU_gf2PolyFree(&Q_temp[i], 0);
	}

	BPU_gf2QcMatrixToSparse(&ctx->code_spec->qcldpc->Q_sparse, &Q_temp_mat, wq);
	BPU_gf2QcMatrixFree(&Q_temp_mat, 0);

	ret = 0;

	BPU_printDebug("calculating Sinv_G = S^-1 x G");
	//Multiplication of matrix S^-1 blocks with matrix G's identity block
	//SG_temp[0,1,2, ..., (n0-2)] contain the main diagonal blocks
	for (i = 0; i < ctx->code_spec->qcldpc->n0 - 1; i++) {
		BPU_gf2PolyCopy(&SG_temp[i], &S_inv_temp[i]);

#if defined(DEBUG_L)
		BPU_printDebug("Sinv_G[%i]: ", i);
		//BPU_printGf2Poly(&SG_temp[i]);
#endif

	}

	//Multiplication of matrix S^-1 blocks with the rest of matrix G
	//SG_temp[(n0-1) ... (n0-1)*2] contain rest of the blocks
	for (i = 0; i < ctx->code_spec->qcldpc->n0 - 1; i++) {
		BPU_gf2PolyMalloc((&SG_temp[i + (ctx->code_spec->qcldpc->n0 - 1)]), ctx->code_spec->qcldpc->m);
		for (j = 0; j < ctx->code_spec->qcldpc->m; j++) {
			if (BPU_gf2PolyGetBit((&S_inv_temp[i]), j))
				BPU_gf2VecXor(&SG_temp[i + (ctx->code_spec->qcldpc->n0 - 1)], &G_temp[i]);
			BPU_gf2PolyMulX(&G_temp[i]);
		}
		BPU_gf2PolyFree(&S_inv_temp[i], 0);
		BPU_gf2PolyFree(&G_temp[i], 0);

#if defined(DEBUG_L)
		BPU_printDebug("Sinv_G[%i]: ", (i + (ctx->code_spec->qcldpc->n0 - 1)));
		//BPU_printGf2Poly(&SG_temp[i]);
#endif
	}

	BPU_printDebug("calculating G_masked = S^-1 x G x Q^-1");
	//G_masked[0,1,2, ..., (n0-2)] contain the main diagonal blocks
	//G_masked[(n0-1) ... (n0-1)*2] contain the rest of the blocks

	//for every row of blocks in SG_temp
	for (i = 0; i < (ctx->code_spec->qcldpc->n0 - 1); i++) {
		BPU_gf2PolyMalloc((&G_masked[i]), ctx->code_spec->qcldpc->m);
		BPU_gf2PolyMalloc((&G_masked[i + (ctx->code_spec->qcldpc->n0 - 1)]), ctx->code_spec->qcldpc->m);
		//for every row in a block
		for (j = 0; j < ctx->code_spec->qcldpc->m; j++) {
			if (BPU_gf2PolyGetBit((&SG_temp[i]), j)) {
				BPU_gf2VecXor(&G_masked[i], &Q_inv_temp[i]);
			}
			if (BPU_gf2PolyGetBit((&SG_temp[i + (ctx->code_spec->qcldpc->n0 - 1)]), j)) {
				BPU_gf2VecXor(&G_masked[i + (ctx->code_spec->qcldpc->n0 - 1)], &Q_inv_temp[ctx->code_spec->qcldpc->n0 - 1]);
			}
			//cyclic shift by 1
			BPU_gf2PolyMulX(&Q_inv_temp[i]);
			BPU_gf2PolyMulX(&Q_inv_temp[ctx->code_spec->qcldpc->n0 - 1]);
		}
		BPU_gf2PolyFree(&Q_inv_temp[i], 0);
		BPU_gf2PolyFree(&SG_temp[i], 0);
		BPU_gf2PolyFree(&SG_temp[i + (ctx->code_spec->qcldpc->n0 - 1)], 0);
	}
	BPU_gf2PolyFree(&Q_inv_temp[ctx->code_spec->qcldpc->n0 - 1], 0);

	BPU_printDebug("copying G_masked into code context");
	BPU_gf2QcMatrixMalloc(&ctx->code_spec->qcldpc->G_masked, (ctx->code_spec->qcldpc->n0 - 1) * 2, ctx->code_spec->qcldpc->m, 0, 0);
	for (i = 0; i < (ctx->code_spec->qcldpc->n0 - 1) * 2; i++) {
		BPU_gf2PolyCopy(&ctx->code_spec->qcldpc->G_masked.matrices[i], &G_masked[i]);
		BPU_gf2PolyFree(&G_masked[i], 0);
	}

	BPU_printDebug("free and exit");
	free(S_temp);
	free(S_inv_temp);
	free(Q_temp);
	free(Q_inv_temp);
	free(SG_temp);
	free(G_masked);
	BPU_gf2PolyFree(&mod, 0);

	return ret;
}

int BPU_QcldpcMaskKeys2(BPU_T_Code_Ctx *ctx, BPU_T_GF2_Poly *G_temp) { 
	BPU_T_GF2_Poly *SG_temp, *G_masked;
	BPU_T_GF2_QC_Matrix Qmat, Qinv, Sinv;
	BPU_T_Perm_Vector *perm;

	SG_temp = (BPU_T_GF2_Poly*)calloc(ctx->code_spec->qcldpc->n0 * (ctx->code_spec->qcldpc->n0-1), sizeof(BPU_T_GF2_Poly));
	G_masked = (BPU_T_GF2_Poly*)calloc(ctx->code_spec->qcldpc->n0 * (ctx->code_spec->qcldpc->n0-1), sizeof(BPU_T_GF2_Poly));

	int32_t i, ret, row, column, block, bit, found;
	//Hamming weight of Matrix Q blocks	
	int32_t wq[] = { 3,4,4,4,3,4,4,4,3 };	
	int32_t *tmp_wq = (int32_t *)malloc(sizeof(wq));

	found = 0;
	row = 0;	
	BPU_printDebug("generating matrices S and S^-1");
	while (!found) {
		ret = BPU_gf2QcMultirowMatrixMalloc(&ctx->code_spec->qcldpc->S, ctx->code_spec->qcldpc->n0-1, ctx->code_spec->qcldpc->n0-1, ctx->code_spec->qcldpc->m);
		if (ret != 0) {
			BPU_printError("Could not allocate memory");
			return 1;
		}
		//generating blocks not on the main diagonal
		for (i = 0; i < (ctx->code_spec->qcldpc->n0 - 1)*(ctx->code_spec->qcldpc->n0 - 1); i++) {
			if (i == row*(ctx->code_spec->qcldpc->n0-1) + row)
				continue;
			//generate polynomial of even Hamming weight, probability of bit being set to 1 is 50%
			BPU_gf2PolyInitRandProb(&ctx->code_spec->qcldpc->S.matrices[i], ctx->code_spec->qcldpc->m, 0.5, 0, 1);
			if (i%(ctx->code_spec->qcldpc->n0-1) == 0)
				row++;
		}
		//generating blocks on the main diagonal
		for (i = 0; i < ctx->code_spec->qcldpc->n0-1; i++)
			//generate polynomial of odd Hamming weight, probability of bit being set to 1 is 50%
			BPU_gf2PolyInitRandProb(&ctx->code_spec->qcldpc->S.matrices[i + i*(ctx->code_spec->qcldpc->n0 - 1)],
				ctx->code_spec->qcldpc->m, 0.5, 0, 0);
		//inverting matrix S
		found = BPU_gf2QcInvMatrix(&Sinv, &ctx->code_spec->qcldpc->S);
		//if no inversion exists, free matrices and start anew
		if (!found) {
			BPU_printWarning("inversion to matrix S not found");
			BPU_gf2QcMatrixFree(&ctx->code_spec->qcldpc->S, 0);			
		}
	}

	found = 0;
	row = 0;
	BPU_printDebug("generating matrices Q and Q^-1");
	while (!found) {
		ret = BPU_gf2QcMultirowMatrixMalloc(&Qmat, ctx->code_spec->qcldpc->n0, ctx->code_spec->qcldpc->n0, ctx->code_spec->qcldpc->m);
		if (ret != 0) {
			BPU_printError("Could not allocate memory");
			return 1;
		}
		//generating blocks not on the main diagonal
		for (i = 0; i < (ctx->code_spec->qcldpc->n0)*(ctx->code_spec->qcldpc->n0); i++) {
			if (i == row*ctx->code_spec->qcldpc->n0 + row)
				continue;			
			BPU_gf2PolyInitRand(&Qmat.matrices[i], ctx->code_spec->qcldpc->m, wq[i], 0);
			if (i%ctx->code_spec->qcldpc->n0 == 0)
				row++;
		}
		//generating blocks on the main diagonal
		for (i = 0; i < ctx->code_spec->qcldpc->n0; i++)			
			BPU_gf2PolyInitRand(&Qmat.matrices[i + i*ctx->code_spec->qcldpc->n0], ctx->code_spec->qcldpc->m, wq[i + i*ctx->code_spec->qcldpc->n0], 0);

		//inverting matrix Q
		found = BPU_gf2QcInvMatrix(&Qinv, &Qmat);
		//if no inversion exists, free matrices and start anew
		if (!found){
			BPU_printWarning("inversion to matrix Q not found");
			BPU_gf2QcMatrixFree(&Qmat, 0);			
			}
	}

	BPU_printDebug("applying permutation to matrix Q"); 
	BPU_permMalloc(&perm, ctx->code_spec->qcldpc->n0);	 
	BPU_permRandomize(perm);  
	BPU_gf2QcPermuteMatrixRows(&Qmat, perm);	 
	BPU_gf2QcPermuteMatrixColumns(&Qinv, perm); 
	
	//permuting rows of the array containing Hamming weight of matrix Q's blocks
	for (i = 0; i < perm->size; i++) {
		for (column = 0; column < Qmat.column_element_count; column++) {
			tmp_wq[i * Qmat.column_element_count + column] =
				wq[perm->elements[i] * Qmat.column_element_count + column];
		}
	}
	for (i = 0; i < Qmat.element_count; i++) {
		wq[i] = tmp_wq[i];		
	}
			
	BPU_permFree(&perm);

	BPU_permMalloc(&perm, ctx->code_spec->qcldpc->n0);	
	BPU_permRandomize(perm);		
	BPU_gf2QcPermuteMatrixColumns(&Qmat, perm);	
	BPU_gf2QcPermuteMatrixRows(&Qinv, perm);
	
	//permuting columns of the array containing Hamming weight of matrix Q's blocks
	for (i = 0; i < perm->size; i++) {
		for (row = 0; row < Qmat.row_element_count; row++) {
			tmp_wq[i + Qmat.row_element_count * row] =
				wq[perm->elements[i] + Qmat.row_element_count * row];
		}
	}
	for (i = 0; i < Qmat.element_count; i++)
		wq[i] = tmp_wq[i];
	
	BPU_permFree(&perm);	
	free(tmp_wq);

	// ***********************************************************************
	// fprintf(stderr, "BPU_QcldpcMaskKeys2\n");
	// fprintf(stderr, " - Qmat.n: %i\n", Qmat.n);
	// fprintf(stderr, " - Qmat.k: %i\n", Qmat.k);
	// fprintf(stderr, " - Qmat.element_count: %i\n\n", Qmat.element_count);

	// BPU_mecsQcldpcTestGHmatrices(&Qmat, &ctx->code_spec->qcldpc->H);

	// BPU_T_GF2_Sparse_Qc_Matrix *H = &ctx->code_spec->qcldpc->H;

	// BPU_T_GF2_Matrix * h;
	// BPU_gf2MatMalloc(&h, H->k, H->n);

	// // int i; 
	// int j;
	// int index;
	// // int bit, column;
	// BPU_T_GF2_Sparse_Poly ROW;
	// for (i = 0; i < H->k; i++) {
	// // for (i = 0; i < 1; i++) {
	// 	BPU_gf2SparseQcMatrixGetRow(&ROW, H, i);

	// 	// fprintf(stderr, "Sparse poly (%i):\n", row.weight);
	// 	for (j = 0; j < ROW.weight; j++) {
	// 		index = ROW.index[j];
	// 		bit = index % 32;	
	// 		column = index / 32;
	// 		h->elements[i][column] |= (1u << (31 - bit));

	// 		// fprintf(stderr, "%4i: %3i - %2i  ", index, column, bit);
	// 		// BPU_printBinaryMsbLn(h->elements[i][column], 32);	
	// 	}
	// 	BPU_gf2SparsePolyFree(&ROW, 0);
	// }

	// BPU_T_GF2_Matrix * q;
	// BPU_gf2MatMalloc(&q, Qmat.n, Qmat.k);

	// for (i = 0; i < Q->k; i++) {
	// 	BPU_gf2SparseQcMatrixGetRow(&row, Q, i);

	// 	for (j = 0; j < row.weight; j++) {
	// 		index = row.index[j];
	// 		bit = i % 32;
	// 		column = i / 32;
	// 		q->elements[index][column] |= (1u << (31 - bit));
	// 	}

	// 	BPU_gf2SparsePolyFree(&row, 0);
	// }

	BPU_printDebug("copying sparse matrix Q into code context");
	BPU_gf2QcMatrixToSparse(&ctx->code_spec->qcldpc->Q_sparse, &Qmat, wq);
	
	// fprintf(stderr, "BPU_QcldpcMaskKeys2\n");
	// fprintf(stderr, " - ctx->code_spec->qcldpc->Q_sparse.n: %i\n", ctx->code_spec->qcldpc->Q_sparse.n);
	// fprintf(stderr, " - ctx->code_spec->qcldpc->Q_sparse.k: %i\n\n", ctx->code_spec->qcldpc->Q_sparse.k);
	// ***********************************************************************


	BPU_printDebug("calculating Sinv_G = S^-1 x G");
	//Multiplication of matrix S^-1 blocks with matrix G's identity block
	for (row = 0; row < ctx->code_spec->qcldpc->n0 - 1; row++) {
		for (column = 0; column < ctx->code_spec->qcldpc->n0 - 1; column++) {
			BPU_gf2PolyCopy(&SG_temp[column + row*ctx->code_spec->qcldpc->n0], &Sinv.matrices[column + row*(ctx->code_spec->qcldpc->n0 - 1)]);

#if defined(DEBUG_L)
			BPU_printDebug("Sinv_G[%i]: ", column + row*ctx->code_spec->qcldpc->n0);
			//BPU_printGf2Poly(&SG_temp[column + row*ctx->code_spec->qcldpc->n0]);
#endif

		}
	}

	//Multiplication of matrix S^-1 blocks with the rest of matrix G
	for (row = 0; row < ctx->code_spec->qcldpc->n0 - 1; row++) {
		BPU_gf2PolyMalloc(&SG_temp[(row + 1)*(ctx->code_spec->qcldpc->n0) - 1], ctx->code_spec->qcldpc->m);
		for (block = 0; block < ctx->code_spec->qcldpc->n0 - 1; block++) {
			for (bit = 0; bit < ctx->code_spec->qcldpc->m; bit++) {
				if (BPU_gf2PolyGetBit((&Sinv.matrices[row * (ctx->code_spec->qcldpc->n0 - 1) + block]), bit))
					BPU_gf2VecXor(&SG_temp[(row + 1)*(ctx->code_spec->qcldpc->n0) - 1], &G_temp[block]);
				BPU_gf2PolyMulX(&G_temp[block]);
			}
		}
	}
	for(i=0;i<ctx->code_spec->qcldpc->n0-1;i++)
        BPU_gf2PolyFree(&G_temp[i],0);
	
	BPU_printDebug("calculating G_masked = S^-1 x G x Q^-1");
	for (row = 0; row < ctx->code_spec->qcldpc->n0 - 1; row++) {
		for (column = 0; column < ctx->code_spec->qcldpc->n0; column++) {
			BPU_gf2PolyMalloc((&G_masked[column + row*ctx->code_spec->qcldpc->n0]), ctx->code_spec->qcldpc->m);
			for (block = 0; block < ctx->code_spec->qcldpc->n0; block++) {
				for (bit = 0; bit < ctx->code_spec->qcldpc->m; bit++) {
					if (BPU_gf2PolyGetBit((&SG_temp[block + row*ctx->code_spec->qcldpc->n0]), bit))
						BPU_gf2VecXor(&G_masked[column + row*ctx->code_spec->qcldpc->n0], (&Qinv.matrices[block * ctx->code_spec->qcldpc->n0 + column]));
					BPU_gf2PolyMulX(&Qinv.matrices[block * ctx->code_spec->qcldpc->n0 + column]);
				}
			}
		}
	}

	BPU_printDebug("copying G_masked into code context");
	BPU_gf2QcMatrixMalloc(&ctx->code_spec->qcldpc->G_masked, ctx->code_spec->qcldpc->n0 * (ctx->code_spec->qcldpc->n0-1), ctx->code_spec->qcldpc->m, 0, 0);
	for (i = 0; i < ctx->code_spec->qcldpc->n0 * (ctx->code_spec->qcldpc->n0-1); i++) {
		BPU_gf2PolyCopy(&ctx->code_spec->qcldpc->G_masked.matrices[i], &G_masked[i]);		
		BPU_gf2PolyFree(&G_masked[i], 0);
		BPU_gf2PolyFree(&SG_temp[i],0);
	}

	// fprintf(stderr, "BPU_QcldpcMaskKeys2\n");
	// fprintf(stderr, " - ctx->code_spec->qcldpc->G_masked.n: %i\n", ctx->code_spec->qcldpc->G_masked.n);
	// fprintf(stderr, " - ctx->code_spec->qcldpc->G_masked.k: %i\n\n", ctx->code_spec->qcldpc->G_masked.k);

	BPU_printDebug("free and exit");
	free(SG_temp);
	free(G_masked);
	BPU_gf2QcMatrixFree(&Qmat, 0);
	BPU_gf2QcMatrixFree(&Qinv, 0);
	BPU_gf2QcMatrixFree(&Sinv, 0);
	return ret;
}

int BPU_mecsQcldpcTestGHmatrices(const BPU_T_GF2_QC_Matrix *G, const BPU_T_GF2_Sparse_Qc_Matrix *H) {
	int element, err = 0;
	BPU_T_GF2_Poly temp;
	BPU_T_GF2_Sparse_Poly row;
	uint32_t i;

	// get I * H[0] + ... + I * H[n0-2]
	BPU_gf2PolyMalloc(&temp, G->n);
	for (i = 0; i < G->element_count; i++) {
		BPU_gf2SparsePolyAdd(&temp, &H->matrices[i]);
	}

	// get other rows
	for (element = 0; element < G->element_count; element++) {
		for (i = 0; i < G->element_size; i++) {
			if (BPU_gf2VecGetBit(&G->matrices[element], i) == 1ul) {
				BPU_gf2SparseQcMatrixGetRow(&row, H, i + G->element_size * (G->element_count));
				BPU_gf2SparsePolyAdd(&temp, &row);
				BPU_gf2SparsePolyFree(&row, 0);
			}
		}
	}
	// check result poly
	for (i = 0; i < temp.elements_in_row; i++) {
		if (temp.elements[i] != 0ul) {
			err++;
			break;
		}
	}

#if defined(DEBUG_L)
	BPU_T_GF2_QC_Matrix GH_result;
	BPU_gf2QcMatrixMalloc(&GH_result, 1, G->element_size, 0, 0);
	BPU_gf2PolyCopy(&GH_result.matrices[0], &temp);
	//    BPU_printGf2QcMatrix(&GH_result);
	BPU_gf2QcMatrixFree(&GH_result, 0);
#endif

	BPU_gf2PolyFree(&temp, 0);

	return err;
}

#endif
