///2203867996232225326
/// This file was automatically generated using the VnV-Matcher executable. 
/// The matcher allows for automatic registration of all VnV plugins and injection 
/// points. Building the matcher requires Clang. If Clang is not available on this machine,
/// Registration code should be written manually. 
/// 

//PACKAGENAME: HYPRE_EX1

#include "VnV.h" 
DECLARESUBPACKAGE(VnVHypre)
const char* getFullRegistrationJson_HYPRE_EX1(){
	 return "{\"Communicator\":{\"docs\":\"\",\"name\":\"mpi\",\"package\":\"VNV\"},\"Conclusion\":\"\",\"Introduction\":\"\",\"SubPackages\":{\"VnVHypre\":{\"docs\":\"\",\"name\":\"VnVHypre\",\"packageName\":\"HYPRE_EX1\"}}}";}

INJECTION_REGISTRATION(HYPRE_EX1){
	REGISTERSUBPACKAGE(HYPRE_EX1,VnVHypre);
	VnV_Declare_Communicator("HYPRE_EX1","VNV","mpi");
	REGISTER_FULL_JSON(HYPRE_EX1, getFullRegistrationJson_HYPRE_EX1);
};



