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

%macro malloc_func 1
	mov  rcx, %1                   ; request %1 bytes
	call malloc                    ; allocate memory
%endmacro

section .data
	floatformat: 		DB 	"%lf",0
	printcheckFormat:	DB	"rel %lf img %lf",4,0
	neg:			DB	-1

SECTION .TEXT
	GLOBAL cmplx_add_s
	GLOBAL cmplx_sub_s
	GlOBAL cmplx_mult_s
	GlOBAL cmplx_div_s

	GlOBAL func_malloc

	extern malloc
	extern free
	extern printf

section .data
	tmp:	DQ 0.0	;global floating point varaiebl
	divisor DQ 0.0  ;global floating point for divide function

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
	mov qword[r8], 0	;*res_real=0
	mov qword[r9], 0	;*res_img=0
	;real_mult




	;(a_img*b_img)
	fld qword [rsi]		;st0=[rsi]=a_img
	fst st1			;st1=st0=a_img
	fld qword [rcx]		;st0=M[rsi]=b_img
	fmul			;st0=a_img*b_img
	fst qword [r8]		;M[r8]=st0=(a_img*b_img)

	;a_real*b_real
	fld qword [rdi]		;st0=M[rdi]=a_real
	fst st1			;st1=st0=a_real
	fld qword [rdx]		;st0=M[rdx]=b_real
	fmul			;st0=a_real*b_real

	;a_real*b_real-(a_img*b_img)
	fst st1			;st1=-(a_img*b_img)
	fld qword [r8]		;st0=a_real*b_real
	fsub			;st0=a_real*b_real+-(a_img*b_img)
	fst qword [r8]		;*res_real=st0=a_real*b_real+-(a_img*b_img)


	;img_mult

	;a_real*b_img
	fld qword [rdi]		;st0=M[rdi]=a_real
	fst st1			;st1=st0=a_real
	fld qword [rcx]		;st0=M[rcx]=b_img
	fmul
	fst qword[r9]		;st2=st0=a_real*b_img

	;a_img*b_real
	fld qword [rsi]		;st0=M[rsi]=a_img
	fst st1			;st1=st0=a_img
	fld qword [rdx]		;st0=M[rdx]=b_real
	fmul

	
	;a_img*b_real+a_real*b_img
	fst st1			;st1=a_img*b_real
	fld qword [r9]		;st0=a_real*b_img
	fadd			;st0=a_a_img*b_real+a_real*b_img
	fst qword [r9]		;*res_img=st0=a_a_img*b_real+a_real*b_img

	funcend

;;important:divide of floating point- st0 always is divisor.
cmplx_div_s:

	funcstart

	;;divisor-(b.real)^2+(b.img)^2
	

	;;b.real^b.real
	fld qword [rdx]		;st0=M[rdx]=b_real
	fst st1			;st1=st0=b_real
	fld qword [rdx]		;st0=M[rdx]=b_real
	fmul			;st0=b_real*b_real
	fst qword [tmp]		;tmp=b_real^2


	;;b.img*b.img
	fld qword [rcx]		;st0=M[rcx]=b_img
	fst st1			;st1=st0=b_img
	fld qword [rcx]		;st0=M[rcx]=b_img
	fmul			;st0=b_img*b_img

	;;(b.real)^2+(b.img)^2
	fst st1			;st1=(b.img)^2
	fld qword [tmp]		;st0=(b.real)^2
	fadd			;st0=(b.real)^2+(b.img)^2
	fst qword [divisor]	;divisor=(b.real)^2+(b.img)^2


	;;real section- a.real*b.real+a.img*b.img

	;(a_img*b_img)
	fld qword [rsi]		;st0=[rsi]=a_img
	fst st1			;st1=st0=a_img
	fld qword [rcx]		;st0=M[rsi]=b_img
	fmul			;st0=a_img*b_img
	fst qword [tmp]		;tmp=st0=(a_img*b_img)

	;a_real*b_real
	fld qword [rdi]		;st0=M[rdi]=a_real
	fst st1			;st1=st0=a_real
	fld qword [rdx]		;st0=M[rdx]=b_real
	fmul			;st0=a_real*b_real

	;a_real*b_real+(a_img*b_img)
	fst st1			;st1=a_real*b_real
	fld qword[tmp]		;st0=tmp=(a_img*b_img)
	fadd			;st0=a_real*b_real+(a_img*b_img)

	;a_real*b_real+(a_img*b_img)/divisor
	fst st1			;st1=a_real*b_real+(a_img*b_img)
	fld qword [divisor]	;st0=divisor=(b.real)^2+(b.img)^2
	fdiv			;st0=a_real*b_real+(a_img*b_img)/divisor
	fst qword [r8]		;res_real=a_real*b_real+(a_img*b_img)/divisor

	;;img section- a.img*b.real-a.real*b.img/divisor

	;a.img*b.real-a.real*b.img

	;a.real*b.img
	fld qword [rdi]		;st0=M[rdi]=a_real
	fst st1			;st1=st0=a_real
	fld qword [rcx]		;st0=M[rcx]=b_img
	fmul			;st0=a_real*b_img
	fst qword[tmp]		;tmp=a.img*b.real

	;a.img*b.real
	fld qword [rsi]		;st0=M[rsi]=a_img
	fst st1			;st1=st0=a_img
	fld qword [rdx]		;st0=M[rdx]=b_real
	fmul			;st0=a.img*b.real

	;a.img*b.real-a.real*b.img
	fst st1			;st1=st0=a_real*b_img
	fld qword [tmp]		;st0=a.img*b.real
	fsub			;st0=a.img*b.real-a_real*b_img

	;a.img*b.real-a.real*b.img/divisor
	fst st1			;st1=a.img*b.real-a_real*b_img
	fld qword [divisor]	;st0=divisor=(b.real)^2+(b.img)^2
	fdiv			;st0=a.img*b.real-a_real*b_img/divisor
	fst qword [r9]		;res_img=a.img*b.real-a_real*b_img/divisor

	funcend

func_malloc:
	funcstart
	malloc_func 8

	mov  qword  [rax],   1      ; write "1" into low 32 bits
	mov  qword  [rax+4], 2      ; write "2" into high 32 bits

	funcend


















