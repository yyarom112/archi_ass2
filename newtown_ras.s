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

%macro iplusplus 0
	push r14
	mov  r14, qword[i]                   	; r14=i
	inc r14                    		
	mov qword[i],r14
	pop r14 
%endmacro

%macro iminusminus 0
	push r14
	mov  r14, qword[i]                   	; r14=i
	dec r14                    		
	mov qword[i],r14 
	pop r14 
%endmacro

%macro init_globals 0
	mov qword [tmp],0
	mov qword [divisor],0
	mov qword [pow],0
	mov qword [i],0
	mov qword [n],0
	mov qword [tmp_real],0
	mov qword [tmp_img],0
%endmacro

%macro callfunc 0
	push rdi
	push rsi
	push rdx
	push rcx
	push rbx
	push r8
	push r9
	push r10
	push r11
	push r12
	push r13
	push r14
	push r15

%endmacro

%macro returnfunc 0
	pop r15
	pop r14
	pop r13
	pop r12
	pop r11
	pop r10
	pop r9
	pop r8
	pop rbx
	pop rcx
	pop rdx
	pop rsi
	pop rdi
%endmacro

section .data


	tmp:	DQ 0.0	;global floating point varaiebl
	divisor: DQ 0.0  ;global floating point for divide function
	pow: DQ 0	;global int for eval_derivative ->pow of argumentS
	i: DQ 0		;global int for eval_derivative ->run index
	n: DQ 0		;global int for sqoort
	tmp_real:DQ 0	;Eveal Poly
	tmp_img: DQ 0	;Eval poly
	;return_val1:DQ 0
	;return_val2:DQ 0
	epsilon: DQ 0
	;order: DQ 0
	index_coeff: DQ 0
	real_coeff: DQ 0			
	img_coeff: DQ 0
	;init_img_coeff: DQ 0	;init point arg- NOT TOUCH-> only in main!!
	;init_real_coeff: DQ 0	;init point arg- NOT TOUCH-> only in main!!
	int_tmp: DQ 0

	;;meseges for scanf & printf
	epsilon_msg:db " epsilon = %lf",0
	order_input:db " order = %d",0
	coeff_input: db " coeff %d = %lf %lf",0
	initial_point_input: db " initial = %lf %lf",0
	print_val: db "value = %d",0
	print_something: db "*************",0
	print_final_result: db "root = %lf %lf",0
	print_tmp: db "tmp = %lf  ",0


section .bss
	return_val1:	resq 1
	return_val2:	resq 1
	order:		resq 1
	div_order:	resq 1
	div_real_arr:	resq 1	
	div_img_arr:	resq 1
	real_arr:	resq 1
	img_arr:	resq 1
	f_real:		resq 1
	f_img:		resq 1
	fd_real:	resq 1
	fd_img:		resq 1
	c_real:		resq 1
	c_img:		resq 1
	init_img_coeff:	resq 1
	init_real_coeff:resq 1

	arr:		resq 1	

SECTION .TEXT
	GLOBAL cmplx_add_s
	GLOBAL cmplx_sub_s
	GlOBAL cmplx_mult_s
	GlOBAL cmplx_div_s
	GlOBAL eval_derivative_s
	GlOBAL func_malloc
	GlOBAL sqroot_s
	GlOBAL make_normal_s
	GlOBAL close_enough_s
	GlOBAL eval_poly_s
	GlOBAL newton_rashford_impl_s
	GLOBAL main



	GlOBAL test_arr_float

	extern malloc
	extern calloc
	extern free
	extern printf
	extern scanf


;argument order: rdi->rsi->rdx->rcx->r8->r9->stack


main:
	;nop
	
	enter 0,0				;enter point
	finit
	funcstart
	;scanf-epslion
	callfunc
	mov rax,1				;num arg for scanf=1
	mov rdi,epsilon_msg			;arg=formaot of get epsilon
	lea rsi, [epsilon]			;return val-> epsilon
	call scanf
	returnfunc

	;scanf-order
	callfunc
	mov rax,1				;num arg for scanf=1
	mov rdi,order_input			;arg=formaot of get order
	lea rsi, [int_tmp]			;return val-> order
	call scanf
	returnfunc
	mov rdi, 0
	mov dword eax,[int_tmp]
	mov dword[order],eax

	;let place to arrays
	mov r14,qword[order]	
	mov rdi,r14			;rdi=order
	dec r14
	add rdi,r14
	add rdi,rdi			;rdi=2order+2(order-1)
	mov rsi,8			;size of one place in arr is sizeof(double)
	call calloc
	mov qword[arr],rax		;arr=return val of calloc

	;;coeff while
	
	mov r9,0
	mov r10, qword[order]		;r10=order
	mov qword[i],r10		;i is num of iterration
	cmp r10,0			;check if the input order is 0
	jl .coeff_while_end		;****maybe need to change to jle
	.coeff_while:
		;;prepare args
		finit
		callfunc
		mov rax,1
		mov rdi,coeff_input	;get the scanf params into rdi
		lea rsi,[index_coeff]
		lea rdx,[real_coeff]
		lea rcx,[img_coeff]
		mov qword[tmp],0
		lea r11,[tmp]
		call scanf
		returnfunc

		;working on the index of insertion ->r9 
		mov r9,qword[order]
		mov r12,qword[index_coeff]
		sub r9,r12
		fld qword[real_coeff]
		fstp qword[arr+8*r9]	;real_coeff is double- size
		
		mov r12, qword[order]
		add r9,r12		; r9= points to the img part of the arr
		fld qword[img_coeff]
		fstp qword[arr+8*r9]

		fldz			;init to zero the args
		fst qword[index_coeff]
		fst qword[real_coeff]		
		fst qword[img_coeff]		

		;prepare next iteration
		iminusminus
		mov r12, qword[i]
		;cmp r12,0		;if i<0 ->end of iterations
		;dec r10
		;cmp r10,0
		jl .coeff_while_end
		jmp .coeff_while

	.coeff_while_end:	

	;get the initial_point_input
	callfunc	
	mov rax,1				;num of args is only 1- 1 line of input left
	mov rdi,initial_point_input		;get the scanf params into rdi
	lea rsi,[init_real_coeff]
	lea rdx,[init_img_coeff]
	call scanf
	returnfunc

	;arrays-pointers define
	mov r10, qword [arr]			;r10=arr
	mov qword[real_arr],r10			;real_arr=*arr=*real_arr
	mov r11, qword [order]			;r11=order
	imul r11,8				;need to jump *8 because double
	add r10,r11
	mov qword[img_arr],r10			;img_arr=*img_arr
	add r10,r11				;after adding to r10 order*2 r10 points to the start of the div_real_arr
	mov qword[div_real_arr], r10
	sub r11,8				;r11=order-1(in double)
	add r10,r11				;r10=arr+2*order+(order-1)
	mov qword[div_img_arr], r10



	;;void newton_rashford_impl(double *real_arr , double* img_arr, double epsilon,int order, double* init_real, double* init_img)
	callfunc				;backup registers
	mov rdi,real_arr
	mov rsi,img_arr
	mov rdx,epsilon
	mov rcx,order
	mov r8, init_real_coeff
	mov r9, init_img_coeff
	call newton_rashford_impl_s
	returnfunc

	;;print the result
	push rdi
	mov rdi,print_final_result
	movsd xmm0,qword[init_real_coeff]
	movsd xmm1,qword[init_img_coeff]
	mov rax, 2
	call printf
	pop rdi
	funcend



	





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
	;;arguments 
	funcstart
	
	mov qword[i],0			;i=0
	mov qword[pow],0		;pow=0

	.for:
		mov r10, qword[i]	;r10=i
		mov r11, r8		;r8=len
		sub r10,r8		;r10=i-len
		cmp r10,0
		jge .enf_for		;if(i-len>=0) jump to end
	
		;pow=len-i
		mov r11, r8		;r8=len
		mov r10, qword[i]	;r12=i
		sub r11,r10		;r11=r11-r10=len-i
		mov qword[pow], r11	;pow=len-i
		
		;;tmp= real_arr[i]*pow;
		mov r10, qword[i]	;r10=i
		mov r11, qword[pow]	;r11=pow
		fld qword[rdi+8*r10]	;st0=*[real_arr+i*sizeof(double))]
		fst st1			;st1=real_arr+i*sizeof(double))
		fild qword[pow]		;st0=pow
		fmul			;st0=*[real_arr+i*sizeof(double))]*pow
		fst qword[tmp]		;tmp=*[real_arr+i*sizeof(double))]*pow
		
		;res_real_arr[i]=tmp
		mov r13,qword[tmp]	;r13=*[real_arr+i*sizeof(double))]*pow
		mov r10, qword[i]	;r10=i
		mov qword[rdx+8*r10],r13;res_real_arr[i]=tmp

		;;tmp= img_arr[i]*pow;
		mov r10, qword[i]	;r10=i
		mov r11, qword[pow]	;r11=pow
		fld qword[rsi+8*r10]	;st0=*[img_arr+i*sizeof(double))]
		fst st1			;st1=img_arr+i*sizeof(double))
		fild qword[pow]		;st0=pow
		fmul			;st0=*[img_arr+i*sizeof(double))]*pow
		fst qword[tmp]		;tmp=*[img_arr+i*sizeof(double))]*pow
		
		;res_img_arr[i]=tmp
		mov r13,qword[tmp]	;r13=*[img_arr+i*sizeof(double))]*pow
		mov r10, qword[i]	;r10=i
		mov qword[rcx+8*r10],r13;res_real_arr[i]=tmp


		;preper the next interation
		iplusplus			;need r14 available
		jmp .for



	.enf_for:	

	funcend




;;for now this func is void and get double num and double* ret_val
sqroot_s:
	;;get one arg inf rdi=n
	funcstart
	
	;;mov qword[n],rdi		;n=rdi=arg1
	;fld qword[n]			;st0=n
	fld qword[rdi]	
	fsqrt				;st0=sqrt(n)
	fst qword[return_val1]			;n=st0=sqrt(n)
	;mov rax, qword[rsi]		;return_val=n=sqrt(n)	
	
	funcend


make_normal_s:
	;;got two arg: rdi=*a_real, rsi=*a_img
	
	funcstart

	;tmp=img*img
	fld qword[rsi]			;st0=*a_img
	fst st1				;st1=*a_img
	fmul				;st0=a_img*a_img
	fld qword[tmp]			;tmp=a_img*a_img
	
	;a_real*a_real
	fld qword[rdi]			;st0=*a_real
	fst st1				;st1=*a_real
	fmul				;st0=a_real*a_real

	;;output=a_real*a_real+a_img*a_img
	fst st1				;st1=a_real*a_real
	fld qword[tmp]			;st0=tmp=a_img*a_img
	fadd				;st0=a_img*a_img+a_real*a_real
	fst qword[tmp]			;tmp=a_img*a_img+a_real*a_real
	
	;;sqrt to tmp
	
	callfunc 
	
	mov rdi, tmp			;arg1=a_img*a_img+a_real*a_real
	mov rsi, n			;output=n
	call sqroot_s
	returnfunc
	;mov qword[return_val1],rax	
	funcend 







close_enough_s:
	;;got 3 arg: rdi=*a_real, rsi=*a_img, rdx=epsilon rcx=return_value

	funcstart
	
	callfunc			;back up the register
	call make_normal_s

	returnfunc			;restore the registers

	fld qword[return_val1]
	fst qword[rcx]


	funcend
	

eval_poly_s:
	;argument list: rdi=*real_arr rsi=*img_arr rdx=z_img rcx=z_img r8=res_real r9=_res_img	

	funcstart

	mov qword[i],1			;i=1
	;*res img=img_arr[0]
	fld qword[rsi]			;st0=img_arr[0]
	fst qword[r9]			;*res_img=img_arr[0]

	;*res_real=real_arr[0]
	fld qword[rdi]			;st0=arr_real[0]
	fst qword[r8]			;*res_real=real_arr[0]
	
	mov qword[tmp_img],0		;tmp_img=0
	mov qword[tmp_real],0		;tmp_real=0
	
	.for:
		;for check cond
		mov r10, qword[i]	;r10=i
		mov r11, qword[pow]	;r11=pow
		sub r10,r11		;r10=i-pow
		cmp r10,0
		jg  .end_for		;if(i-pow>0) jump to end

		; cmplx_mult(*res_real,*res_img,*z_real,*z_img,&tmp_real,&tmp_img)
		callfunc			;back up registers
		
		;;prepare args
		mov rdi,r8			;arg1=rdi=res_real
		mov rsi,r9			;arg2=rsi=res_img
		mov r8,	tmp_real		;arg5=tmp_real
		mov r9, tmp_img			;arg6=tmp_img

		call cmplx_mult_s
		
		returnfunc			;restore the register
		
		;tmp=real_arr[i]
		mov r10, qword[i]		;r10=i
		mov r11,qword[rdi+8*r10]	;r11=real_arr[i]
		mov qword[tmp], r11		;tmp=real_arr[i]

		;n=img_arr[i]
		mov r10, qword[i]		;r10=i
		mov r11,qword[rsi+8*r10]	;r11=img_arr[i]
		mov qword[n], r11		;n=img_arr[i]

		;cmplx_add(tmp_real,tmp_img,real_arr[i],img_arr[i],res_real,res_img);
		
		callfunc			;back up registers
		
		mov rdi,tmp_real		;arg1=rdi=tmp_real
		mov rsi,tmp_img			;arg2=rsi=tmp_img
		mov rdx,tmp			;arg3=rdx=real_arr[i]
		mov rcx,n			;arg4=rcx=img_arr[i]
		
		call cmplx_add_s
		
		returnfunc

		mov qword[tmp_real],0		;tmp_real=0
		mov qword[tmp_img],0		;tmp_img=0

		iplusplus			;i++
		jmp .for


	.end_for:	
		funcend	



newton_rashford_impl_s:
	;list: rdi=*real_arr rsi=*img_arr rdi=*epsilon rcx=*order r8=*init_real r9=*init_img
	
	funcstart
	

	;;eval_derivative(real_arr,img_arr,div_real_arr,div_img_arr,order);
	callfunc				;backup register
	;real_arr and img_arr in place
	mov rdx,div_real_arr
	mov rcx,div_img_arr
	mov r8, qword[order]
	call eval_derivative_s
	returnfunc
	
	
	;;double f_real,f_img,fd_real,fd_img,c_real,c_img;
	;initial values	
	mov r10,0
	mov qword[f_real],r10
	mov qword[f_img],r10
	mov qword[fd_real],r10
	mov qword[fd_img],r10
	mov qword[c_real],r10
	mov qword[c_img],r10

	;; eval_poly(real_arr,img_arr,init_real,init_img,&f_real,&f_img,order);
	callfunc				;backup register
	;real_arr and img_arr are in place
	mov rdx,init_real_coeff
	mov rcx,init_img_coeff
	mov r8,f_real
	mov r9,f_img
	mov r10,qword[order]			;because pow changes from order to order-1 in different cases
	mov qword[pow],r10
	call eval_poly_s
	returnfunc

	;;eval_poly(div_real_arr,div_img_arr,init_real,init_img,&fd_real,&fd_img,order-1);
	callfunc				;backup register
	mov rdi,div_real_arr
	mov rsi,div_img_arr
	mov rdx,init_real_coeff
	mov rcx,init_img_coeff
	mov r8,fd_real
	mov r9,fd_img
	mov r10,qword[order]			;because pow changes from order to order-1 in different cases
	dec r10
	mov qword[div_order],r10
	mov qword[pow],r10
	call eval_poly_s
	returnfunc

	;;cmplx_div(f_real,f_img,fd_real,fd_img,&c_real,&c_img);
	callfunc				;backup register
	mov rdi,f_real
	mov rsi,f_img
	mov rdx,fd_real
	mov rcx,fd_img
	mov r8,c_real
	mov r9,c_img
	call cmplx_div_s
	returnfunc

	.while:
		;check cond
		callfunc				;backup register
		mov rdi,c_real
		mov rsi,c_img
		call make_normal_s
		mov qword[tmp],rax			;move the result of normal in tmp
		returnfunc
		fld qword[epsilon]
		fst st1					;st1=st0=epsilon
		fld qword[return_val1]				;st0=tmp=normal
		fsub					;st0=normal-epsilon
		fst qword[tmp]				;tmp=normal-epsilon
		fldz
		fst st1
		fld qword[tmp]
		fcom
		jl .while_end

		;;the body of while
		;; eval_poly(real_arr,img_arr,init_real,init_img,&f_real,&f_img,order);
		callfunc				;backup register
		;real_arr and img_arr are in place
		mov rdx,init_real_coeff
		mov rcx,init_img_coeff
		mov r8,f_real
		mov r9,f_img
		mov r10,qword[order]			;because pow changes from order to order-1 in different cases
		mov qword[pow],r10
		call eval_poly_s
		returnfunc

		;;eval_poly(div_real_arr,div_img_arr,init_real,init_img,&fd_real,&fd_img,order-1);
		callfunc				;backup register
		mov rdi,div_real_arr
		mov rsi,div_img_arr
		mov rdx,init_real_coeff
		mov rcx,init_img_coeff
		mov r8,fd_real
		mov r9,fd_img
		mov r10,qword[order]			;because pow changes from order to order-1 in different cases
		dec r10
		mov qword[div_order],r10
		mov qword[pow],r10
		call eval_poly_s
		returnfunc

		;;cmplx_div(f_real,f_img,fd_real,fd_img,&c_real,&c_img);
		callfunc				;backup register
		mov rdi,f_real
		mov rsi,f_img
		mov rdx,fd_real
		mov rcx,fd_img
		mov r8,c_real
		mov r9,c_img
		call cmplx_div_s
		returnfunc

		;; *init_real-=c_real;
		fld qword[c_real]
		fst st1
		fld qword[init_real_coeff]
		fsub
		fst qword[init_real_coeff]

		;; *init_img-=c_img;
		fld qword[c_img]
		fst st1
		fld qword[init_img_coeff]
		fsub
		fst qword[init_img_coeff]
		
		jmp .while


	.while_end:
		funcend




	 

;***********************************




;----------------------------------------------------------------------------------------------------------------------------------------

				;st0=a_img*a_img
test_arr_float:
	funcstart
	mov r9,0
	mov r11,rsi
	.for:
		;;cmp r11,0   ;if for terminate - jump to end of func
		;;je .end_for
		;;mov [rdi+r9*8],rdx	;res_img_arr[i]=50.2

		;;prepare next iteration
		dec r11			;num_of_iters--
		inc r9			;i++
		jmp .for	
