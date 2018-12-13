/*
 This file is part of BitPunch
 Copyright (C) 2014-2015 Frantisek Uhrecky <frantisek.uhrecky[what here]gmail.com>

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

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>

#include <bitpunch/crypto/hash/sha512.h>
#include <bitpunch/code/qcldpc/qcldpc.h>
#include <bitpunch/code/qcmdpc/qcmdpc.h>

void vypis(BPU_T_Mecs_Ctx * ctx);

int main(int argc, char **argv) {
	srand(time(NULL));

	BPU_T_Mecs_Ctx *ctx = NULL;
	BPU_T_UN_Mecs_Params params;
	
	// p, n0, hwt(H), hwt(e)
	// m, n0, w, t
	BPU_mecsInitParamsQcldpc(&params, 8192, 3, 39, 34);
	BPU_mecsInitCtx(&ctx, &params, BPU_EN_MECS_BASIC_QCLDPC2);
	BPU_mecsGenKeyPair(ctx);

	// BPU_T_GF2_Vector *ct, *pt_in, *pt_out;
	// BPU_gf2VecMalloc(&ct, ctx->ct_len);
	// BPU_gf2VecMalloc(&pt_in, ctx->pt_len);
	// BPU_gf2VecMalloc(&pt_out, ctx->pt_len);

	// int numExp;
	// int sposob1 = 0;
	// int sposob2 = 0;
	// int pocet_experimentov = 10;
	// for (numExp = 0; numExp < pocet_experimentov; numExp++) {
	// 	BPU_gf2VecRand(pt_in, 0);
	// 	BPU_mecsEncrypt(ct, pt_in, ctx);

	// 	BPU_mecsDecrypt(pt_out, ct, ctx);
	// 	// BPU_mecsQcmdpcDecrypt(pt_out, ct, ctx);

	// 	if (!BPU_gf2VecCmp(pt_in, pt_out))
	// 		sposob1++;
	// }

	// fprintf(stderr, "Pocet experimentov: %i\n\n", pocet_experimentov);
	// fprintf(stderr, "Pocet uspechov QC-LDPC sposobu: %i\n", sposob1);
	// fprintf(stderr, "Uspesnost: %.2f %% \n", (double) sposob1 / pocet_experimentov * 100);
	// fprintf(stderr, "Pocet uspechov QC-MDPC sposobu: %i\n", sposob2);
	// fprintf(stderr, "Uspesnost: %.2f %% \n", (double) sposob2 / pocet_experimentov * 100);



	// ***********************************************************************
	BPU_T_GF2_Sparse_Qc_Matrix * H = &(ctx->code_ctx->code_spec->qcldpc->H);
	BPU_T_GF2_Sparse_Qc_Matrix * Q = &(ctx->code_ctx->code_spec->qcldpc->Q_sparse);
	BPU_T_GF2_QC_Matrix * G_masked = &(ctx->code_ctx->code_spec->qcldpc->G_masked);

	// BPU_printGf2SparseQcMatrix(H);
	// BPU_printGf2SparseQcMatrix(Q);
	// ***********************************************************************


	// ***********************************************************************
	BPU_T_GF2_Matrix * h;
	BPU_gf2MatMalloc(&h, H->k, H->n);

	int i, j, k;
	int index, bit, column;
	BPU_T_GF2_Sparse_Poly row;
	for (i = 0; i < H->k; i++) {
	// for (i = 0; i < 1; i++) {
		BPU_gf2SparseQcMatrixGetRow(&row, H, i);

		// fprintf(stderr, "Sparse poly (%i):\n", row.weight);
		for (j = 0; j < row.weight; j++) {
			index = row.index[j];
			bit = index % 32;	
			column = index / 32;
			// h->elements[i][column] |= (1u << (31 - bit));
			h->elements[i][column] |= (1u << (bit));

			// fprintf(stderr, "%4i: %3i - %2i  ", index, column, bit);
			// BPU_printBinaryMsbLn(h->elements[i][column], 32);	
		}
		BPU_gf2SparsePolyFree(&row, 0);
	}

	// fprintf(stderr, "Matrix size: %dx%d\n", h->k, h->n);
	// for (j = 0; j <= h->elements_in_row - 1; j++) {
	// 	BPU_printBinaryMsbLn(h->elements[0][j], h->element_bit_size);
	// }

	// BPU_gf2SparseQcMatrixGetRow(&row, H, 0);
	// BPU_printGf2SparsePoly(&row);
	// BPU_gf2SparsePolyFree(&row, 0);
	// ***********************************************************************


	// ***********************************************************************
	// BPU_T_GF2_Matrix * q;
	// BPU_gf2MatMalloc(&q, Q->n, Q->k);

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

	// for (i = 0; i < q->k; i++) {
	// 	fprintf(stderr, "%i-", i);
	// 	fprintf(stderr, "%i ", q->elements[i][0] & (1u << 31));
	// }

	// BPU_gf2SparseQcMatrixGetRow(&row, Q, 0);
	// BPU_printGf2SparsePoly(&row);
	// BPU_gf2SparsePolyFree(&row, 0);
	// ***********************************************************************


	// ***********************************************************************
	BPU_T_GF2_Matrix * q;
	BPU_gf2MatMalloc(&q, Q->n * 3, Q->n * 3);

	// fprintf(stderr, "q->k: %i\n", q->k);
	// fprintf(stderr, "q->n: %i\n", q->n);
	// fprintf(stderr, "q->elements_in_row: %i\n", q->elements_in_row);
	// fprintf(stderr, "q->element_bit_size: %i\n\n", q->element_bit_size);

	// for (i = 0; i < 2; i++) {
	for (i = 0; i < Q->k; i++) {
		BPU_gf2SparseQcMatrixGetRow(&row, Q, i);
		// fprintf(stderr, "Sparse poly (%i):\n", row.weight);

		for (j = 0; j < row.weight; j++) {
			index = row.index[j];
			bit = index % 32;
			column = index / 32;
			// q->elements[(i % Q->n) + (i / q->k) * Q->n][column + (q->elements_in_row / 3) * ((i / Q->element_size) % 3)] |= (1u << (31 - bit));
			q->elements[(i % Q->n) + (i / q->k) * Q->n][column + (q->elements_in_row / 3) * ((i / Q->element_size) % 3)] |= (1u << (bit));
		
			// fprintf(stderr, "%4i: %3i - %2i  ", index, column, bit);
			// BPU_printBinaryMsbLn(q->elements[(i % Q->n) + (i / q->k) * Q->n][column + (q->elements_in_row / 3) * ((i / Q->element_size) % 3)], 32);	
		}
		
		BPU_gf2SparsePolyFree(&row, 0);
	}

	// fprintf(stderr, "Matrix size: %dx%d\n", q->k, q->n);
	// for (j = 0; j <= q->elements_in_row - 1; j++) {
	// 	BPU_printBinaryMsbLn(q->elements[0][j], q->element_bit_size);
	// }

	// BPU_gf2SparseQcMatrixGetRow(&row, Q, 0);
	// fprintf(stderr, "%i. riadok: ", 0);
	// BPU_printGf2SparsePoly(&row);
	// BPU_gf2SparsePolyFree(&row, 0);

	// BPU_gf2SparseQcMatrixGetRow(&row, Q, Q->n);
	// fprintf(stderr, "%i. riadok: ", Q->n);
	// BPU_printGf2SparsePoly(&row);
	// BPU_gf2SparsePolyFree(&row, 0);

	// BPU_gf2SparseQcMatrixGetRow(&row, Q, Q->n*2);
	// fprintf(stderr, "%i. riadok: ", Q->n*2);
	// BPU_printGf2SparsePoly(&row);
	// BPU_gf2SparsePolyFree(&row, 0);
	// ***********************************************************************


	// ***********************************************************************
	// REALNE SME SPRAVILI Q*H^T !!!

	BPU_T_GF2_Matrix * HQ;
	BPU_gf2MatMalloc(&HQ, q->k, h->n);
	BPU_gf2MulMat(&HQ, q, h);

	fprintf(stderr, "HQ->k: %i\n", HQ->k);
	fprintf(stderr, "HQ->n: %i\n\n", HQ->n);
	// fprintf(stderr, "HQ->elements_in_row: %i\n", HQ->elements_in_row);
	// fprintf(stderr, "HQ->element_bit_size: %i\n\n", HQ->element_bit_size);

	// BPU_printGf2Mat(HQ);

	// BPU_gf2SparseQcMatrixGetRow(&row, H, 0);
	// BPU_printGf2SparsePoly(&row);
	// BPU_gf2SparsePolyFree(&row, 0);

	// fprintf(stderr, "Matrix size: %dx%d\n", HQ->k, HQ->n);
	// for (j = 0; j <= HQ->elements_in_row - 1; j++) {
	// 	BPU_printBinaryMsbLn(HQ->elements[0][j], HQ->element_bit_size);
	// }
	// ***********************************************************************


	// ***********************************************************************
	BPU_T_GF2_Matrix * g_masked;
	BPU_gf2MatMalloc(&g_masked, G_masked->k * 2, G_masked->k * 3);

	fprintf(stderr, "g_masked->k: %i\n", g_masked->k);
	fprintf(stderr, "g_masked->n: %i\n\n", g_masked->n);
	// fprintf(stderr, "g_masked->elements_in_row: %i\n", g_masked->elements_in_row);
	// fprintf(stderr, "g_masked->element_bit_size: %i\n", g_masked->element_bit_size);

	// vypis(ctx);	
	// BPU_printGf2Poly(&G_masked->matrices[0]);
	// BPU_printGf2PolyForMatrix(&G_masked->matrices[0]);
	// BPU_printGf2QcMatrix(G_masked);

	// fprintf(stderr, "G_masked->matrices[0].elements[0]: %i\n", G_masked->matrices[0].elements[0]);
	// BPU_printBinaryLsbLn(G_masked->matrices[0].elements[0], 32);
	// BPU_printBinaryMsbLn(G_masked->matrices[0].elements[0], 32);
	
	int l = 0;
	BPU_T_GF2_Poly temp;

	for (i = 0; i < 2; i++) {
		for (j = 0; j < 3; j++) {
			fprintf(stderr, "i: %i, j: %i\n", i, j);
			
			BPU_gf2PolyCopy(&temp, &G_masked->matrices[i*3 + j]);			
			for (k = 0; k < G_masked->k; k++) {
				for (l = 0; l < g_masked->elements_in_row / 3; l++) {
					g_masked->elements[k + (i * G_masked->k)][j * (g_masked->elements_in_row / 3) + l] = temp.elements[l];
				}

				BPU_gf2PolyMulX(&temp);
			}
			BPU_gf2PolyFree(&temp, 0);
		} 
	}

	// BPU_printGf2Mat(g_masked);

	// BPU_printGf2Poly(&G_masked->matrices[0]);
	// for (j = 0; j < (g_masked->elements_in_row / 3) - 250; j++) {
	// 	BPU_printBinaryLsb(g_masked->elements[j][0], g_masked->element_bit_size);
	// 	fprintf(stderr, " ");
	// 	BPU_printBinaryLsbLn(g_masked->elements[j][0 + 1], g_masked->element_bit_size);
	// }
	
	// fprintf(stderr, "\n");	

	// BPU_printGf2Poly(&G_masked->matrices[3]);
	// for (j = 0; j < (g_masked->elements_in_row / 3); j++) {
	// 	BPU_printBinaryMsb(g_masked->elements[G_masked->k][j], g_masked->element_bit_size);
	// }
	// ***********************************************************************


	// ***********************************************************************
	fprintf(stderr, "MO:\n\n");	

	BPU_T_GF2_Matrix * MO;
	BPU_gf2MatMalloc(&MO, g_masked->k, HQ->n);
	BPU_gf2MulMat(&MO, g_masked, HQ);

	fprintf(stderr, "MO->k: %i\n", MO->k);
	fprintf(stderr, "MO->n: %i\n\n", MO->n);

	// for (j = 0; j < MO->elements_in_row; j++) {
	// 	BPU_printBinaryMsbLn(MO->elements[0][j], MO->element_bit_size);
	// }

	BPU_printGf2Mat(MO);

	// ***********************************************************************
	
	
	// BPU_gf2VecFree(&ct);
	// BPU_gf2VecFree(&pt_in);
	// BPU_gf2VecFree(&pt_out);
	
	BPU_mecsFreeCtx(&ctx);
	BPU_mecsFreeParamsQcldpc(&params);

	return 0;
}

void vypis(BPU_T_Mecs_Ctx * ctx) {
	fprintf(stderr, "ctx...\n");

	// fprintf(stderr, "typ: %i\n", ctx->type);
	fprintf(stderr, "pt_len: %i\n", ctx->pt_len);
	fprintf(stderr, "ct_len: %i\n\n", ctx->ct_len);

	fprintf(stderr, "t: %i\n", ctx->code_ctx->t);
	// fprintf(stderr, "e: %p\n", ctx->code_ctx->e);
	// fprintf(stderr, "type: %i\n", ctx->code_ctx->type);
	fprintf(stderr, "msg_len: %i\n", ctx->code_ctx->msg_len);
	fprintf(stderr, "code_len: %i\n\n", ctx->code_ctx->code_len);

	// fprintf(stderr, "math_ctx: %i\n", ctx->code_ctx->math_ctx);
	// fprintf(stderr, "math_ctx->ord: %i\n", ctx->code_ctx->math_ctx->ord);
	// fprintf(stderr, "math_ctx->mod: %i\n", ctx->code_ctx->math_ctx->mod);	
	// fprintf(stderr, "math_ctx->mod_deg: %i\n", ctx->code_ctx->math_ctx->mod_deg);

	fprintf(stderr, "m: %i\n", ctx->code_ctx->code_spec->qcldpc->m);
	fprintf(stderr, "n0: %i\n", ctx->code_ctx->code_spec->qcldpc->n0);
	fprintf(stderr, "w: %i\n\n", ctx->code_ctx->code_spec->qcldpc->w);

	fprintf(stderr, "H.n: %i\n", ctx->code_ctx->code_spec->qcldpc->H.n);
	fprintf(stderr, "H.k: %i\n", ctx->code_ctx->code_spec->qcldpc->H.k);
	fprintf(stderr, "H.isVertical: %i\n", ctx->code_ctx->code_spec->qcldpc->H.isVertical);
	fprintf(stderr, "H.element_size: %i\n", ctx->code_ctx->code_spec->qcldpc->H.element_size);
	fprintf(stderr, "H.element_count: %i\n\n", ctx->code_ctx->code_spec->qcldpc->H.element_count);
	// fprintf(stderr, "H.matrices[0].weight: %i\n", ctx->code_ctx->code_spec->qcldpc->H.matrices[0].weight);
	// fprintf(stderr, "H.matrices[0].index[0]: %i\n\n", ctx->code_ctx->code_spec->qcldpc->H.matrices[0].index[0]);

	fprintf(stderr, "Q_sparse.n: %i\n", ctx->code_ctx->code_spec->qcldpc->Q_sparse.n);
	fprintf(stderr, "Q_sparse.k: %i\n", ctx->code_ctx->code_spec->qcldpc->Q_sparse.k);
	fprintf(stderr, "Q_sparse.isVertical: %i\n", ctx->code_ctx->code_spec->qcldpc->Q_sparse.isVertical);
	fprintf(stderr, "Q_sparse.element_size: %i\n", ctx->code_ctx->code_spec->qcldpc->Q_sparse.element_size);
	fprintf(stderr, "Q_sparse.element_count: %i\n\n", ctx->code_ctx->code_spec->qcldpc->Q_sparse.element_count);	
	// fprintf(stderr, "Q_sparse.matrices[0].weight: %i\n", ctx->code_ctx->code_spec->qcldpc->Q_sparse.matrices[0].weight);
	// fprintf(stderr, "Q_sparse.matrices[0].index[0]: %i\n\n", ctx->code_ctx->code_spec->qcldpc->Q_sparse.matrices[0].index[0]);

	// fprintf(stderr, "S.n: %i\n", ctx->code_ctx->code_spec->qcldpc->S.n);
	// fprintf(stderr, "S.k: %i\n", ctx->code_ctx->code_spec->qcldpc->S.k);
	// fprintf(stderr, "S.isVertical: %i\n", ctx->code_ctx->code_spec->qcldpc->S.isVertical);
	// fprintf(stderr, "S.element_size: %i\n", ctx->code_ctx->code_spec->qcldpc->S.element_size);
	// fprintf(stderr, "S.element_count: %i\n", ctx->code_ctx->code_spec->qcldpc->S.element_count);
	// fprintf(stderr, "S.is_I_appended: %i\n", ctx->code_ctx->code_spec->qcldpc->S.is_I_appended);
	// fprintf(stderr, "S.row_element_count: %i\n", ctx->code_ctx->code_spec->qcldpc->S.row_element_count);
	// fprintf(stderr, "S.column_element_count: %i\n", ctx->code_ctx->code_spec->qcldpc->S.column_element_count);
	// fprintf(stderr, "S.matrices->len: %i\n\n", ctx->code_ctx->code_spec->qcldpc->S.matrices->len);
	// fprintf(stderr, "S.matrices->elements: %i\n", ctx->code_ctx->code_spec->qcldpc->S.matrices->elements);
	// fprintf(stderr, "S.matrices->elements_in_row: %i\n", ctx->code_ctx->code_spec->qcldpc->S.matrices->elements_in_row);
	// fprintf(stderr, "S.matrices->element_bit_size: %i\n\n", ctx->code_ctx->code_spec->qcldpc->S.matrices->element_bit_size);
	
	fprintf(stderr, "G_masked.n: %i\n", ctx->code_ctx->code_spec->qcldpc->G_masked.n);
	fprintf(stderr, "G_masked.k: %i\n", ctx->code_ctx->code_spec->qcldpc->G_masked.k);
	fprintf(stderr, "G_masked.isVertical: %i\n", ctx->code_ctx->code_spec->qcldpc->G_masked.isVertical);
	fprintf(stderr, "G_masked.element_size: %i\n", ctx->code_ctx->code_spec->qcldpc->G_masked.element_size);
	fprintf(stderr, "G_masked.element_count: %i\n", ctx->code_ctx->code_spec->qcldpc->G_masked.element_count);
	fprintf(stderr, "G_masked.is_I_appended: %i\n", ctx->code_ctx->code_spec->qcldpc->G_masked.is_I_appended);
	fprintf(stderr, "G_masked.row_element_count: %i\n", ctx->code_ctx->code_spec->qcldpc->G_masked.row_element_count);
	fprintf(stderr, "G_masked.column_element_count: %i\n", ctx->code_ctx->code_spec->qcldpc->G_masked.column_element_count);
	fprintf(stderr, "G_masked.matrices->len: %i\n", ctx->code_ctx->code_spec->qcldpc->G_masked.matrices->len);
	// fprintf(stderr, "G_masked.matrices->elements: %i\n", ctx->code_ctx->code_spec->qcldpc->G_masked.matrices->elements);
	fprintf(stderr, "G_masked.matrices->elements_in_row: %i\n", ctx->code_ctx->code_spec->qcldpc->G_masked.matrices->elements_in_row);
	fprintf(stderr, "G_masked.matrices->element_bit_size: %i\n\n", ctx->code_ctx->code_spec->qcldpc->G_masked.matrices->element_bit_size);
}
