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
#include "qcldpctypes.h"

void BPU_qcldpcFreeSpec(BPU_T_Qcldpc_Spec *spec) {
	
	BPU_gf2QcMatrixFree(&spec->G_masked, 0);
	BPU_gf2SparseQcMatrixFree(&spec->H, 0);
	BPU_gf2SparseQcMatrixFree(&spec->Q_sparse, 0);	
	BPU_gf2QcMatrixFree(&spec->S, 0);
		
	free(spec);
}

int BPU_qcldpcInitParams(BPU_T_Qcldpc_Params **params, const uint16_t m, const uint16_t n0, const uint16_t w, const uint16_t t) {
	*params = (BPU_T_Qcldpc_Params*) calloc(sizeof(BPU_T_Qcldpc_Params), 1);

	if (!params) {
		BPU_printError("Can't init Code params");

		return -1;
	}
	(*params)->m = m;
	(*params)->n0 = n0;
	(*params)->w = w;
	(*params)->t = t;
	
	return 0;
}

void BPU_qcldpcFreeParams(BPU_T_Qcldpc_Params **params) {
	if (*params) {
		free(*params);
	}
	*params = NULL;
}
