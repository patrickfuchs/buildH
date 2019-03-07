function normalize(vec) result(nvec)!{{{
real(kind=rk),intent(in) :: vec(:)
real(kind=rk) :: nvec(1:3)
nvec=vec/sqrt(sum(vec**2))
end function normalize!}}}
    
function v2q(v,theta) result(q)
!Create a quaternion from the vector v and angle theta
real(kind=rk),intent(in) :: v(:),theta ! Rad
real(kind=rk) :: q(4)
q=[cos(theta/2),sin(theta/2)*normalize(v)]
end function v2q

function RV(q) result(m)
!Rotational matrix: B.Stevensson et. al. 2011
real(kind=rk),intent(in) :: q(:) !quat
real(kind=rk) :: m(3,3) !quat
!First col
m(1,1)=q(1)**2+q(2)**2-q(3)**2-q(4)**2
m(2,1)=2*(q(2)*q(3)+q(1)*q(4))
m(3,1)=2*(q(2)*q(4)-q(1)*q(3))
!Second col
m(1,2)=2*(q(2)*q(3)-q(1)*q(4))
m(2,2)=q(1)**2-q(2)**2+q(3)**2-q(4)**2
m(3,2)=2*(q(3)*q(4)+q(1)*q(2))
!Third col
m(1,3)=2*(q(2)*q(4)+q(1)*q(3))
m(2,3)=2*(q(3)*q(4)-q(1)*q(2))
m(3,3)=q(1)**2-q(2)**2-q(3)**2+q(4)**2
end function RV


subroutine add_hydrogen(helpers,hcoor,imol,func,bondlength,theta)
    integer(kind=ik),intent(in) :: helpers(:),func,imol
    integer(kind=ik) :: i
    real(kind=rk) ::&
    v1(3),v2(3),v3(3),v4(3),v5(3),u(3),bondlength,thetal
    real(kind=rk),intent(inout) :: hcoor(:)
    real(kind=rk),optional :: theta
    if(.not.present(theta) .and. func==3)stop 'CH2x needs rotation angle'
    select case(func)
    case(1) ! CH
        v1=atom(helpers(1))%coor(:,imol)
        v2=0._rk
        do i=2,size(helpers)
        ! Center of the nomalized vectors
        v2=v2+normalize(atom(helpers(i))%coor(:,imol)-v1)
        end do
        ! Normalize by number and shift it back to v1
        v2=v2/(size(helpers)-1)+v1
        hcoor=bondlength*normalize(v1-v2)+v1
    case(2) !CH double bond
        v1=atom(helpers(1))%coor(:,imol)
        v2=normalize(atom(helpers(2))%coor(:,imol)-v1)
        v3=normalize(atom(helpers(3))%coor(:,imol)-v1)
        ! thetal is the angle 2pi - C-C-C devided by 2
        ! to ensure equal (C-C-H) angles from both directions
        thetal=pi*(2-valenceangle(helpers(1),&
                helpers(2),&
                helpers(3),imol)/180._rk)/2
        u=normalize(cross_product(v2,v3))
        ! RV returns a rotational 3x3 matrix from a quaternion
        ! v2q turns a vector into a quaternion with angle theta
        ! Rotate by matrix multiplication RV*V3 where V3 is rotated
        hcoor=bondlength*normalize(matmul(RV(v2q(u,thetal)),v3))+v1
    case(3) !CH2
        thetal=theta*pi/180._rk
        v1=atom(helpers(1))%coor(:,imol)
        v2=normalize(atom(helpers(2))%coor(:,imol)-v1)
        v3=normalize(atom(helpers(3))%coor(:,imol)-v1)
        ! Perpendicular to the C-C-C plane
        v4=normalize(cross_product(v3,v2))
        ! Vector to rotate around. In the plane.
        u=normalize(v2-v3) 
        ! The vector to be rotated theta/2, perpendicular to u and v4
        v4=normalize(cross_product(v4,u))
        ! Rotate the new v4 around u by thetal 
        hcoor=bondlength*normalize(matmul(RV(v2q(u,thetal)),v4))+v1
    case(4) !CH3e
        thetal=acos(-1._rk/3) ! 109.47 deg in radians
        v1=atom(helpers(1))%coor(:,imol)
        v2=normalize(atom(helpers(2))%coor(:,imol)-v1)
        v3=normalize(atom(helpers(3))%coor(:,imol)-v1)
        u=normalize(cross_product(v3,v2))
        hcoor=bondlength*normalize(matmul(RV(v2q(u,thetal)),v2))+v1
    case(5) !CH3r
        thetal=120._rk*pi/180._rk
        v1=atom(helpers(1))%coor(:,imol)
        v2=atom(helpers(2))%coor(:,imol)
        v3=atom(helpers(3))%coor(:,imol) !First hydrogen added by CH3e
        u=normalize(v2-v1) ! C-C bond to rotate around
        v4=normalize(v3-v1) 
        hcoor=bondlength*normalize(matmul(RV(v2q(u,thetal)),v4))+v1
    case(6) !CH3s
        thetal=-120._rk*pi/180._rk
        v1=atom(helpers(1))%coor(:,imol)
        v2=atom(helpers(2))%coor(:,imol)
        v3=atom(helpers(3))%coor(:,imol) !First hydrogen added by CH3e
        u=normalize(v2-v1) ! C-C bond to rotate around
        v4=normalize(v3-v1)
        hcoor=bondlength*normalize(matmul(RV(v2q(u,thetal)),v4))+v1
    end select
end subroutine add_hydrogen