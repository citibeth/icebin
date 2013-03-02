namespace glint2 {
namespace modele {


// Example from GEOM_B.f:
// 
// print *,'GEOM_SPECS', i_0h, i_1h, j_0h, j_1h,i_0, i_1, j_0, j_1, j_0s, j_1s
// Num MPI Processes:4
//  GEOM_SPECS     1    72     0    12     1    72     1    11     1
//  GEOM_SPECS     1    72    35    47     1    72    36    46     1
//  GEOM_SPECS     1    72    11    24     1    72    12    23     1
//  GEOM_SPECS     1    72    23    36     1    72    24    35     1
// 
// // From model/MPI_Support/dd2d_utils.f
// 
//       integer :: i_0h,i_1h,j_0h,j_1h,i_0,i_1,j_0,j_1,j_0s,j_1s
// 
//       i_0h = grid%i_strt_halo
//       i_1h = grid%i_stop_halo
//       j_0h = grid%j_strt_halo
//       j_1h = grid%j_stop_halo
//       i_0 = grid%i_strt
//       i_1 = grid%i_stop
//       j_0 = grid%j_strt
//       j_1 = grid%j_stop
//       j_0s = grid%j_strt_skp
//       j_1s = grid%j_stop_skp
// 
// 
// From model/MPI_Support/dd2d_utils.f
// 
//          ! Parameters for Global domain
//         INTEGER :: IM_WORLD     ! Number of Longitudes
//         INTEGER :: JM_WORLD     ! Number of latitudes
//          ! Parameters for local domain
//         INTEGER :: I_STRT       ! Begin local domain longitude index
//         INTEGER :: I_STOP       ! End   local domain longitude index
//         INTEGER :: J_STRT       ! Begin local domain latitude  index
//         INTEGER :: J_STOP       ! End   local domain latitude  index
//         INTEGER :: J_STRT_SKP   ! Begin local domain exclusive of S pole
//         INTEGER :: J_STOP_SKP   ! End   local domain exclusive of N pole
//         INTEGER :: ni_loc       ! for transpose
//          ! Parameters for halo of local domain
//         INTEGER :: I_STRT_HALO  ! Begin halo longitude index
//         INTEGER :: I_STOP_HALO  ! End   halo longitude index
//         INTEGER :: J_STRT_HALO  ! Begin halo latitude  index
//         INTEGER :: J_STOP_HALO  ! End   halo latitude  index


// ------------------------------------------------
class ModelEIndex : public glint2::LocalIndex {
public:
	int i_f;		// Fortran-style, starting with 1
	int j_f;		// Fortran-style, starting with 1

	ModelEIndex(int _i_f, int _j_f) : i_f(_i_f), j_f(_j_f) {}
};

// ------------------------------------------------
class ModelEDomain : public glint2::GridDomain {

	int im_f, jm_f;
	int i0h_f, i1h_f, j0h_f, j1h_f;
	int i0_f, i1_f, j0_f, j1_f;
	int j0s_f, j1s_f;

	ModelEDomain(
		// Info about the global grid
		int _im_f, int _jm_f,

		// Info about the local grid (C-style indices)
		int _i0h_f, int _i1h_f, int _j0h_f, int _j1h_f,
		int _i0_f, int _i1_f, int _j0_f, int _j1_f,
		int _j0s_f, int _j1s_f) : GridDomain(2),

		im_f(_im_f), jm_f(_jm_f),
		i0h_f(_i0h_f), i1h_f(_i1h_f), j0h_f(_j0h_f), j1h_f(_j1h_f),
		i0_f(_i0_f), i1_f(_i1_f), j0_f(_j0_f), j1_f(_j1_f),
		j0s_f(_j0s_f), j1s_f(_j1s_f) {}


	/** Given a global index (C-style 0...ndata()-1), returns a local
	index for this MPI node.  Returns false if this is outside our halo. */
	bool global_to_local(int gindex_c, int *lindex)
	{
		// Decompose index into i & j (zero-based)
		int j_c = gindex / im;
		int i_c = gindex - im * j_c;

		// Convert to 1-based indexing
		lindex[0] = i_c + 1;
		lindex[1] = j_c + 1;

		return true;
	}

	bool in_domain(int *lindex)
	{
	}

	bool in_domain(int *lindex)
	{
	}


	bool in_domain(int gindex_c)
	{
		int j_c = gindex / im;
		int i_c = gindex - im * j_c;
		int j_f = j_c + 1;
		int i_f = i_c + 1;

		return (j_f >= j0_f && j_f <= j1_f);
	}

	bool in_halo(int gindex_c)
	{
		int j_c = gindex / im;
		int i_c = gindex - im * j_c;
		int j_f = j_c + 1;
		int i_f = i_c + 1;

		return (j_f >= j0h_f && j_f <= j1h_f);
	}
};
// ------------------------------------------------





}}
