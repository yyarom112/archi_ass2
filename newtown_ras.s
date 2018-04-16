%macro funcstart 0
	push		rbp
	mov		rbp, rsp
	finit
%endmacro

%macro funcend 0
	mov		rsp, rbp
	pop		rbp
	ret
%endmacro
section .data
	floatformat: 		DB 	"%lf",0
	printcheckFormat:	DB	"rel %lf img %lf",4,0
	neg:			DB	-1

SECTION .TEXT
	GLOBAL cmplx_add_s
	GLOBAL cmplx_sub_s
	GlOBAL cmplx_mult_s
	extern malloc
	extern free
	extern printf

;argument order: rdi->rsi->rdx->rcx->r8->r9->stack


cmplx_add_s:
	funcstart
	mov qword[r8], 0	;*res_real=0
	mov qword[r9], 0	;*res_img=0

	;real_add
	fld qword [rdi]		;st0=M[rdi]=a_real
	fst st1			;st1=st0=a_real
	fld qword [rdx]		;st0=M[rdx]=b_real
	fadd			;st0+=st1
	fst qword [r8]		;*res_real=st0=a_real+b_real

	;img_add
	fld qword [rsi]		;st0=[rsi]=a_img
	fst st1			;st1=st0=a_img
	fld qword [rcx]		;st0=M[rsi]=b_img
	fadd			;st0+=st1
	fst qword [r9]		;*res_img=st0=a_img+b_img

	funcend

cmplx_sub_s:
	funcstart
	mov qword[r8], 0	;*res_real=0
	mov qword[r9], 0	;*res_img=0

	;real_sub
	fld qword [rdi]		;st0=M[rdi]=a_real
	fst st1			;st1=st0=a_real
	fld qword [rdx]		;st0=M[rdx]=b_real
	fsub			;st0-=st1
	fst qword [r8]		;*res_real=st0=a_real-b_real

	;img_sub
	fld qword [rsi]		;st0=[rsi]=a_img
	fst st1			;st1=st0=a_img
	fld qword [rcx]		;st0=M[rsi]=b_img
	fsub			;st0-=st1
	fst qword [r9]		;*res_img=st0=a_img-b_img

	funcend

cmplx_mult_s:
	funcstart
	
	;real_mult

	;a_real*b_real
	fld qword [rdi]		;st0=M[rdi]=a_real
	fst st1			;st1=st0=a_real
	fld qword [rdx]		;st0=M[rdx]=b_real
	fmul

	fst st2			;st2=st0=a_real*b_real

	;-(a_img*b_img)
	fld qword [rsi]		;st0=[rsi]=a_img
	fst st1			;st1=st0=a_img
	fld qword [rcx]		;st0=M[rsi]=b_img
	fmul			;st0=a_img*b_img

	fst st1			;st1=a_img*b_img
	fld qword [neg]		;st0=-1
	fmul			;st0=-(a_img*b_img)


	;a_real*b_real+-(a_img*b_img)
	fst st1			;st1=-(a_img*b_img)
	fld st2			;st0=a_real*b_real
	fadd			;st0=a_real*b_real+-(a_img*b_img)
	fst qword [r8]		;*res_real=st0=a_real*b_real+-(a_img*b_img)


	;img_mult

	;a_real*b_img
	fld qword [rdi]		;st0=M[rdi]=a_real
	fst st1			;st1=st0=a_real
	fld qword [rcx]		;st0=M[rcx]=b_img
	fmul
	fst st2			;st2=st0=a_real*b_img

	;a_img*b_real
	fld qword [rsi]		;st0=M[rsi]=a_img
	fst st1			;st1=st0=a_img
	fld qword [rdx]		;st0=M[rdx]=b_real
	fmul

	
	;a_img*b_real+a_real*b_img
	fst st1			;st1=a_img*b_real
	fld st2			;st0=a_real*b_img
	fadd			;st0=a_a_img*b_real+a_real*b_img
	fst qword [r9]		;*res_img=st0=a_a_img*b_real+a_real*b_img

	funcend



