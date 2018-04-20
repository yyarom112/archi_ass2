%macro funcstart 0
	push		rbp
	mov		rbp, rsp
	finit
	;init_globals
%endmacro

%macro funcend 0
	mov		rsp, rbp
	pop		rbp
	ret
%endmacro

%macro malloc_macro 1
	mov  rcx, %1                   ; request %1 bytes
	call malloc                    ; allocate memory
%endmacro

%macro init_globals 0
	mov qword [tmp],0
	mov qword [divisor],0
	mov qword [pow],0
	mov qword [i],0
%endmacro

section .data
	floatformat: 		DB 	"%lf",0
	printcheckFormat:	DB	"rel %lf img %lf",4,0
	neg:			DB	-1
	tmp:	DQ 0.0	;global floating point varaiebl
	divisor: DQ 0.0  ;global floating point for divide function
	pow: DQ 0	;global int for eval_derivative ->pow of argumentS
	i: DQ 0		;global int for eval_derivative ->run index
	

SECTION .TEXT
	GLOBAL cmplx_add_s
	GLOBAL cmplx_sub_s
	GlOBAL cmplx_mult_s
	GlOBAL cmplx_div_s
	GlOBAL eval_derivative_s
	GlOBAL func_malloc


	GlOBAL test_arr_float

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
	fsub			;st0=a_real*b_real-(a_img*b_img)
	fst qword [r8]		;*res_real=st0=a_real*b_real-(a_img*b_img)

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
	malloc_macro 8

	mov  qword  [rax],   1      ; write "1" into low 32 bits
	mov  qword  [rax+4], 2      ; write "2" into high 32 bits

	funcend


eval_derivative_s:
	;; r9=i , r10=pow , r11=number of iterations , r12=tmp_register
	;; example for torder:  x^2+3x+4
	;; assumption - len= the highest order in the funs, not the size of array
	
	funcstart
;; r9=i , r10=pow , r11=number of iterations , r12=tmp_register
	;; example for torder:  x^2+3x+4
	;; assumption - len= the highest order in the funs, not the size of array
	
	funcstart
	;mov r9,0	;i=0
	mov r10,0	;pow=0	
	mov r11,r8	;r11=r8=len (for the iteration)
	.for:
		cmp r11,0   ;if for terminate - jump to end of func
		je .end_for
		
		;;calc pow
		mov r10,r8		;r10=pow= len
		;sub r10,r9		;r10=pow=len-i
		;;TEST
		mov qword[tmp],r10
		mov r9 ,[tmp]
		funcend


		;;calc real part
		mov r12,[rdi+r9*8] 	;r12=real_arr[i]
		mov [tmp],r12		;tmp=r12=real_arr[i]
		mov [pow],r10		;pow=r10=len-i
		fld qword [tmp]		;st0=global[tmp]=real_arr[i]
		fst st1			;st1=st0=real_arr[i]
		fld qword [pow]		;st0=global[pow]=pow
		fmul			;st0=real_arr[i]*pow
		fst qword[tmp]		;tmp=real_arr[i]*pow
		mov r12,qword [tmp]	;r12=tmp=real_arr[i]*pow
		mov [rdx+r9*8],r12	;res_real_arr[i]=real_arr[i]*pow



		;;calc img part
		mov r12,0
		mov r12,[rsi+r9*8]	;r12=img_arr[i]
		mov [tmp],r12		;tmp=r12=img_arr[i]
		mov [pow],r10		;pow=r10=len-i
		fld qword [tmp]		;st0=global[tmp]=img_arr[i]
		fst st1			;st1=st0=img_arr[i]
		fld qword [pow]		;st0=global[pow]=pow
		fmul			;st0=img_arr[i]*pow
		fst qword[tmp]		;tmp=img_arr[i]*pow
		mov r12,qword [tmp]	;r12=tmp=img_arr[i]*pow

		mov [rcx+r9*8],r12	;res_img_arr[i]=img_arr[i]*pow

		;;prepare next iteration
		dec r11			;num_of_iters--
		inc r9			;i++
		jmp .for
	.end_for:
	funcend








test_arr_float:
	funcstart
	mov r9,0
	mov r11,rsi
	.for:
		cmp r11,0   ;if for terminate - jump to end of func
		je .end_for
		mov [rdi+r9*8],rdx	;res_img_arr[i]=50.2

		;;prepare next iteration
		dec r11			;num_of_iters--
		inc r9			;i++
		jmp .for		



	.end_for:
		funcend

		;;test
		;;mov rdx, 3
		;mov rcx, r12
		;;funcend




