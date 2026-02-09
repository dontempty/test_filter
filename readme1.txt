LES : dynamic model
x,y non-uniform : MG
x uniform : CS
inflow routine
iteration routine
make_ini rouinte 
Blasius routine 

tbl_absoft : make 라는 명령을 사용하여 컴파일
	리눅스 기계
	IMSL 루틴
	랜덤함수

	cosine transform or multigrid
	1)include POISSON_CS_IMSL.F or POISSON_MG_IMSL.F
	2)READUP 에서 해당하는 부분 활성화

	참고사항
	02_main
		PARAM.H
			MT=0	: 위상평균 아님
			MT .ne. 0 	: 위상평균


	iteration : 	GETUP에서 CALL UHCALC_I 를 호출
	decoupling:                      CALL UHCALC 를 호출

	LES : tbl.set -> ILES=1
		EDDY_VISCOS_DM : dynamic model : Germano et al., Lilly
			test filtering : Simpson rule &  Delta_tilde/Delta_bar=2
		EDDY_VISCOS_SM : Smagorinsky model

	@ 연구노트 참조 : CODE NOTE : LES + MG

tbl_compaq
	컴팩기계 or  비주얼 포트란(포아송 부분에서 include .90를 활성화)
	CXML 루틴
	랜덤함수

	POISSON_CS_CXML.F
	POISSON_MG_CXML.F