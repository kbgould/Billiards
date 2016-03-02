struct VECTOR {
	double x;
	double y;
	double z;
};

// UTILITY VECTORS

VECTOR* UP_AXIS = new VECTOR;
VECTOR* Tmp_Vec_1 = new VECTOR;
VECTOR* Tmp_Vec_2 = new VECTOR;

////////////////////
//// TABLE WALLS ///
////////////////////

VECTOR* BL_to_BR1 = new VECTOR;
VECTOR* BL_to_BR2 = new VECTOR;

VECTOR* BR_to_MR1 = new VECTOR;
VECTOR* BR_to_MR2 = new VECTOR;

VECTOR* MR_to_TR1 = new VECTOR;
VECTOR* MR_to_TR2 = new VECTOR;

VECTOR* TR_to_TL1 = new VECTOR;
VECTOR* TR_to_TL2 = new VECTOR;

VECTOR* TL_to_ML1 = new VECTOR;
VECTOR* TL_to_ML2 = new VECTOR;

VECTOR* ML_to_BL1 = new VECTOR;
VECTOR* ML_to_BL2 = new VECTOR;

////////////////////////
///// POCKET WALLS /////
////////////////////////

// BOTTOM LEFT
VECTOR* BL_Lower1 = new VECTOR;
VECTOR* BL_Lower2 = new VECTOR;

VECTOR* BL_Upper1 = new VECTOR;
VECTOR* BL_Upper2 = new VECTOR;
	
// BOTTOM RIGHT CORNER
VECTOR* BR_Lower1 = new VECTOR;
VECTOR* BR_Lower2 = new VECTOR;

VECTOR* BR_Upper1 = new VECTOR;
VECTOR* BR_Upper2 = new VECTOR;

// TOP RIGHT CORNER
VECTOR* TR_Lower1 = new VECTOR;
VECTOR* TR_Lower2 = new VECTOR;

VECTOR* TR_Upper1 = new VECTOR;
VECTOR* TR_Upper2 = new VECTOR;

// TOP LEFT CORNER
VECTOR* TL_Lower1 = new VECTOR;
VECTOR* TL_Lower2 = new VECTOR;

VECTOR* TL_Upper1 = new VECTOR;
VECTOR* TL_Upper2 = new VECTOR;

// MID LEFT POCKET
VECTOR* ML_Lower1 = new VECTOR;
VECTOR* ML_Lower2 = new VECTOR;

VECTOR* ML_Upper1 = new VECTOR;
VECTOR* ML_Upper2 = new VECTOR;

// MID RIGHT POCKET
VECTOR* MR_Lower1 = new VECTOR;
VECTOR* MR_Lower2 = new VECTOR;

VECTOR* MR_Upper1 = new VECTOR;
VECTOR* MR_Upper2 = new VECTOR;

/// LARGE WALLS
VECTOR* TopInt = new VECTOR;
VECTOR* RightInt = new VECTOR;
VECTOR* LeftInt = new VECTOR;
VECTOR* BottomInt = new VECTOR;
VECTOR* BL_Low_Int = new VECTOR;
VECTOR* BL_High_Int = new VECTOR;
VECTOR* BR_Low_Int = new VECTOR;
VECTOR* BR_High_Int = new VECTOR;
VECTOR* MR_Low_Int = new VECTOR;
VECTOR* MR_High_Int = new VECTOR;
VECTOR* ML_Low_Int = new VECTOR;
VECTOR* ML_High_Int = new VECTOR;
VECTOR* TL_Low_Int = new VECTOR;
VECTOR* TL_High_Int = new VECTOR;
VECTOR* TR_Low_Int = new VECTOR;
VECTOR* TR_High_Int = new VECTOR;

VECTOR* mid = new VECTOR;
VECTOR* mid_to_B1 = new VECTOR;
VECTOR* mid_to_B2 = new VECTOR;
VECTOR* vc = new VECTOR;
VECTOR* vcn = new VECTOR;
VECTOR* u_to_B1 = new VECTOR;
VECTOR* u_to_B2 = new VECTOR;
VECTOR* prj_v1vc = new VECTOR;
VECTOR* prj_v1vcn = new VECTOR;
VECTOR* prj_v2vc = new VECTOR;
VECTOR* prj_v2vcn = new VECTOR;

VECTOR* thisPocket = new VECTOR;

