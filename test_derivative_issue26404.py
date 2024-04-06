#Generation :
import sympy as sm
from sympy import (cos, sin, Matrix, symbols, zeros)
import sympy.physics.mechanics as me
from random import randrange
from sympy.physics.mechanics import (dynamicsymbols)

def kane_with_derivatives(n=1, cart_force=True, joint_torques=False):
    r"""
        Took kane function but changed some types

    """
    if n <= 0:
        raise ValueError('The number of links must be a positive integer.')
    m = list(sm.symbols('m:{}'.format(n + 1)))
    l = list(sm.symbols('l:{}'.format(n)))
    
    q = [sm.Function('q'+str(i)) (m[i]) for i in range(n+1)]
    u = [q[i].diff((m[i])) for i in range(n+1)]

    g, t = sm.symbols('g t')

    I = me.ReferenceFrame('I')
    O = me.Point('O')
    O.set_vel(I, 0)

    P0 = me.Point('P0')
    P0.set_pos(O, q[0] * I.x)
    P0.set_vel(I, u[0] * I.x)
    Pa0 = me.Particle('Pa0', P0, m[0])

    frames = [I]
    points = [P0]
    particles = [Pa0]
    forces = [(P0, -m[0] * g * I.y)]
    kindiffs = [q[0].diff(t) - u[0]]

    possible = [m, l, q, u]
    matrix=[]
    for i in range(n) :
        new_line=[]
        for j in range(n) :
            s=0
            for k in range(4):
                i=randrange(4)
                s+=possible[i][randrange(200)%len(possible[i])]*possible[(i+1)%4][randrange(200)%len(possible[(i+1)%4])]
            new_line+=[s]
        matrix+=[new_line]
        
    return matrix



##############################################################################################################
#kane 20
##############################################################################################################

def avg_timing(func, l0, m0, l1, m1, l2, m2, l3, m3, l4, m4, l5, m5, l6, m6, l7, m7,\
    l8, m8, l9, m9, l10, m10, l11, m11, l12, m12, l13, m13, l14,\
    m14, l15, m15, l16, m16, l17, m17, l18, m18, l19, m19, m20, g, q0, u0, q1,\
    u1, q2, u2, q3, u3, q4, u4, q5, u5, q6, u6, q7, u7, q8, u8,\
    q9, u9, q10, u10, q11, u11, q12, u12, q13, u13, q14, u14, q15, u15,\
    q16, u16, q17, u17, q18, u18, q19, u19, q20, u20,F, T1):
    import time
    m = 0
    for i in range(100):
        temps = time.time()
        _ = func(l0, m0, l1, m1, l2, m2, l3, m3, l4, m4, l5, m5, l6, m6, l7, m7,\
        l8, m8, l9, m9, l10, m10, l11, m11, l12, m12, l13, m13, l14,\
        m14, l15, m15, l16, m16, l17, m17, l18, m18, l19, m19, m20,g,q0, u0,\
                 q1, u1, q2, u2, q3, u3, q4, u4, q5, u5, q6, u6, q7, u7, q8, u8,\
        q9, u9, q10, u10, q11, u11, q12, u12, q13, u13, q14, u14, q15, u15,\
        q16, u16, q17, u17, q18, u18, q19, u19, q20, u20,F, T1)
        m += time.time()-temps
    print(m/100)


l0, m0, l1, m1, l2, m2, l3, m3, l4, m4, l5, m5, l6, m6, l7, m7,\
    l8, m8, l9, m9, l10, m10, l11, m11, l12, m12, l13, m13, l14,\
    m14, l15, m15, l16, m16, l17, m17, l18, m18, l19, m19, m20 =\
    sm.symbols("l0 m0 l1 m1 l2 m2 l3 m3 l4 m4 l5 m5 l6 m6 l7 m7 l8 m8 l9 m9 l10 m10 l11 m11 l12 m12 l13 m13 l14 m14 l15 m15 l16 m16 l17 m17 l18 m18 l19 m19 m20")

g = sm.symbols("g")

#q symbols
q0, q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, q14, q15, q16, q17, q18, q19, q20 =\
sm.Function('q0')(m0), sm.Function('q1')(m1), sm.Function('q2')(m2), sm.Function('q3')(m3),\
sm.Function('q4')(m4), sm.Function('q5')(m5), sm.Function('q6')(m6), sm.Function('q7')(m7),\
sm.Function('q8')(m8), sm.Function('q9')(m9), sm.Function('q10')(m10), sm.Function('q11')(m11),\
sm.Function('q12')(m12), sm.Function('q13')(m13), sm.Function('q14')(m14), sm.Function('q15')(m15),\
sm.Function('q16')(m16), sm.Function('q17')(m17), sm.Function('q18')(m18), sm.Function('q19')(m19),\
sm.Function('q20')(m20)

#u symbols
u0, u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15, u16, u17, u18, u19, u20=\
q0.diff(m0), q1.diff(m1), q2.diff(m2), q3.diff(m3), q4.diff(m4), q5.diff(m5), q6.diff(m6),\
q7.diff(m7), q8.diff(m8), q9.diff(m9), q10.diff(m10), q11.diff(m11), q12.diff(m12),\
q13.diff(m13), q14.diff(m14), q15.diff(m15), q16.diff(m16), q17.diff(m17), q18.diff(m18),\
q19.diff(m19), q20.diff(m20)

F, T1 = dynamicsymbols("F T1")

#generate matrix


kane1 = kane_with_derivatives(20)

avg_timing(sm.lambdify((l0, m0, l1, m1, l2, m2, l3, m3, l4, m4, l5, m5, l6, m6, l7, m7,\
    l8, m8, l9, m9, l10, m10, l11, m11, l12, m12, l13, m13, l14,\
    m14, l15, m15, l16, m16, l17, m17, l18, m18, l19, m19, m20,\
                 g, q0, u0, q1, u1, q2, u2, q3, u3, q4, u4, q5, u5, q6, u6, q7, u7, q8, u8,\
    q9, u9, q10, u10, q11, u11, q12, u12, q13, u13, q14, u14, q15, u15,\
    q16, u16, q17, u17, q18, u18, q19, u19, q20, u20, F, T1 ), kane1, cse=True),\
       0, 0 ,1, 1, 2, 9, 0, 1, 2 , 0, 1, 2, 13, 1,\
       45,7,8,9,3,1,2,0,4,5,6,7,8,9,4,5,6,1,2,3,0,4,5,78,9,6,2,3,0,2,1,5,4,8,9,\
       7,8,9,4,5,6,2,1,3,4,5,6,4,8,9,7,6,5,9,8,7,8,9,4,5,6,4,2,1,3,1,5,2,3,0,5,2)


##############################################################################################################
# kane 40
##############################################################################################################
def avg_timing(func, l0, m0, l1, m1, l2, m2, l3, m3, l4, m4, l5, m5, l6, m6, l7, m7, l8, m8, l9, m9,\
    l10, m10, l11, m11, l12, m12, l13, m13, l14, m14, l15, m15, l16, m16 ,l17, m17,\
    l18, m18, l19, m19, l20 ,m20 ,l21, m21, l22, m22, l23, m23, l24, m24, l25, m25,\
    l26, m26, l27, m27 ,l28, m28, l29, m29, l30, m30, l31, m31, l32, m32, l33, m33,\
    l34, m34, l35, m35, l36 ,m36, l37, m37 ,l38 ,m38, l39, m39,  m40, g, q0, u0, q1, u1,\
               q2, u2, q3, u3, q4, u4, q5, u5, q6, u6, q7, u7, q8, u8, q9, u9, q10, u10,\
    q11, u11, q12, u12, q13, u13, q14, u14, q15, u15, q16, u16, q17, u17, q18, u18, q19,\
    u19, q20, u20, q21, u21, q22, u22, q23, u23, q24, u24, q25, u25, q26, u26, q27, u27, q28,\
    u28, q29, u29, q30, u30, q31, u31, q32, u32, q33, u33, q34, u34, q35, u35, q36, u36, q37, \
    u37, q38, u38, q39, u39, q40, u40, F, T1):
    import time
    m = 0
    for i in range(100):
        temps = time.time()
        _ = func(l0, m0, l1, m1, l2, m2, l3, m3, l4, m4, l5, m5, l6, m6, l7, m7, l8, m8, l9, m9,\
    l10, m10, l11, m11, l12, m12, l13, m13, l14, m14, l15, m15, l16, m16 ,l17, m17,\
    l18, m18, l19, m19, l20 ,m20 ,l21, m21, l22, m22, l23, m23, l24, m24, l25, m25,\
    l26, m26, l27, m27 ,l28, m28, l29, m29, l30, m30, l31, m31, l32, m32, l33, m33,\
    l34, m34, l35, m35, l36 ,m36, l37, m37 ,l38 ,m38, l39, m39,  m40, g, q0, u0, q1, u1,\
    q2, u2, q3, u3, q4, u4, q5, u5, q6, u6, q7, u7, q8, u8, q9, u9, q10, u10,\
    q11, u11, q12, u12, q13, u13, q14, u14, q15, u15, q16, u16, q17, u17, q18, u18, q19,\
    u19, q20, u20, q21, u21, q22, u22, q23, u23, q24, u24, q25, u25, q26, u26, q27, u27, q28,\
    u28, q29, u29, q30, u30, q31, u31, q32, u32, q33, u33, q34, u34, q35, u35, q36, u36, q37, \
    u37, q38, u38, q39, u39, q40, u40, F, T1)
        m += time.time() - temps
    print(m/100)


l0, m0, l1, m1, l2, m2, l3, m3, l4, m4, l5, m5, l6, m6, l7, m7, l8, m8, l9, m9,\
    l10, m10, l11, m11, l12, m12, l13, m13, l14, m14, l15, m15, l16, m16 ,l17, m17,\
    l18, m18, l19, m19, l20 ,m20 ,l21, m21, l22, m22, l23, m23, l24, m24, l25, m25,\
    l26, m26, l27, m27 ,l28, m28, l29, m29, l30, m30, l31, m31, l32, m32, l33, m33,\
    l34, m34, l35, m35, l36 ,m36, l37, m37 ,l38 ,m38, l39, m39,  m40 =\
    sm.symbols("l0 m0 l1 m1 l2 m2 l3 m3 l4 m4 l5 m5 l6 m6 l7 m7 l8 m8 l9 m9 l10 m10 l11 m11 l12 m12 l13 m13 l14 m14 l15 m15 l16 m16 l17 m17 l18 m18 l19 m19 l20 m20 l21 m21 l22 m22 l23 m23 l24 m24 l25 m25 l26 m26 l27 m27 l28 m28 l29 m29 l30 m30 l31 m31 l32 m32 l33 m33 l34 m34 l35 m35 l36 m36 l37 m37 l38 m38 l39 m39 m40")

g = sm.symbols("g")

q0, q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, q14, q15, q16, q17, q18, q19, q20,\
    q21, q22, q23, q24, q25, q26, q27, q28, q29, q30, q31, q32, q33, q34, q35, q36, q37, q38,\
    q39, q40=\
sm.Function('q0')(m0), sm.Function('q1')(m1), sm.Function('q2')(m2), sm.Function('q3')(m3),\
                       sm.Function('q4')(m4), sm.Function('q5')(m5), sm.Function('q6')(m6),\
sm.Function('q7')(m7), sm.Function('q8')(m8), sm.Function('q9')(m9), sm.Function('q10')(m10),\
sm.Function('q11')(m11), sm.Function('q12')(m12), sm.Function('q13')(m13), sm.Function('q14')(m14),\
sm.Function('q15')(m15), sm.Function('q16')(m16), sm.Function('q17')(m17), sm.Function('q18')(m18),\
sm.Function('q19')(m19), sm.Function('q20')(m20), sm.Function('q21')(m21), sm.Function('q22')(m22),\
sm.Function('q23')(m23), sm.Function('q24')(m24), sm.Function('q25')(m25), sm.Function('q26')(m26),\
sm.Function('q27')(m27), sm.Function('q28')(m28), sm.Function('q29')(m29), sm.Function('q30')(m30),\
sm.Function('q31')(m31), sm.Function('q32')(m32), sm.Function('q33')(m33), sm.Function('q34')(m34),\
sm.Function('q35')(m35), sm.Function('q36')(m36), sm.Function('q37')(m37), sm.Function('q38')(m38),\
sm.Function('q39')(m39), sm.Function('q40')(m40)

u0, u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15, u16, u17, u18, u19, u20,\
    u21, u22, u23, u24, u25, u26, u27, u28, u29, u30, u31, u32, u33, u34, u35, u36, u37, u38,\
    u39, u40=\
q0.diff(m0), q1.diff(m1), q2.diff(m2), q3.diff(m3), q4.diff(m4), q5.diff(m5), q6.diff(m6),\
q7.diff(m7), q8.diff(m8), q9.diff(m9), q10.diff(m10), q11.diff(m11), q12.diff(m12),\
q13.diff(m13), q14.diff(m14), q15.diff(m15), q16.diff(m16), q17.diff(m17), q18.diff(m18),\
q19.diff(m19), q20.diff(m20), q21.diff(m21), q22.diff(m22), q23.diff(m23), q24.diff(m24), \
q25.diff(m25), q26.diff(m26), q27.diff(m27), q28.diff(m28), q29.diff(m29), q30.diff(m30), \
q31.diff(m31), q32.diff(m32), q33.diff(m33), q34.diff(m34), q35.diff(m35), q36.diff(m36),\
q37.diff(m37), q38.diff(m38), q39.diff(m39), q40.diff(m40)

F, T1 = dynamicsymbols("F T1")

#generate matrix


kane1 = kane_with_derivatives(40)

avg_timing(sm.lambdify((l0, m0, l1, m1, l2, m2, l3, m3, l4, m4, l5, m5, l6, m6,\
    l7, m7, l8, m8, l9, m9,\
    l10, m10, l11, m11, l12, m12, l13, m13, l14, m14, l15, m15, l16, m16 ,l17, m17,\
    l18, m18, l19, m19, l20 ,m20 ,l21, m21, l22, m22, l23, m23, l24, m24, l25, m25,\
    l26, m26, l27, m27 ,l28, m28, l29, m29, l30, m30, l31, m31, l32, m32, l33, m33,\
    l34, m34, l35, m35, l36 ,m36, l37, m37 ,l38 ,m38, l39, m39,  m40,q0, u0, q1, u1, q2,\
    u2, q3, u3, q4, u4, q5, u5, q6, u6, q7, u7, q8, u8, q9, u9, q10, u10,\
    q11, u11, q12, u12, q13, u13, q14, u14, q15, u15, q16, u16, q17, u17, q18, u18, q19,\
u19, q20, u20, q21, u21, q22, u22, q23, u23, q24, u24, q25, u25, q26, u26, q27, u27, q28,\
u28, q29, u29, q30, u30, q31, u31, q32, u32, q33, u33, q34, u34, q35, u35, q36, u36, q37, \
u37, q38, u38, q39, u39, q40, u40 ,g, F, T1  ), kane1, cse=True),\
       0, 0 ,1, 1, 2, 9, 0, 1, 2 , 0, 1, 2, 13, 1,\
       45,7,8,9,3,1,2,0,4,5,6,7,8,9,4,5,6,1,2,3,0,4,5,78,9,6,2,3,0,2,1,5,4,8,9,\
       7,8,9,4,5,6,2,1,3,4,5,6,4,8,9,7,6,5,9,8,7,8,9,4,5,6,4,2,1,3,1,5,2,3,0,5,\
       2, 0 ,1, 1, 2, 9, 0, 1, 2 , 0, 1, 2, 13, 1,\
       45,7,8,9,3,1,2,0,4,5,6,7,8,9,4,5,6,1,2,3,0,4,5,78,9,6,2,3,0,2,1,5,4,8,9,\
       7,8,9,4,5,6,2,1,3,4,5,6,4,8,9,7,6,5,9,8,7,8,9,4,5,6,4,2,1,3,1,5)




##############################################################################################################
# kane 100
##############################################################################################################


def avg_timing(func, l0, m0, l1, m1, l2, m2, l3, m3, l4, m4, l5, m5, l6, m6, l7, m7, l8, m8, l9,\
    m9, l10, m10, l11, m11,\
    l12, m12, l13, m13, l14, m14, l15, m15, l16, m16, l17, m17, l18, m18, l19, m19, l20, m20, l21,\
    m21, l22, m22, l23, m23, l24, m24, l25, m25, l26, m26, l27, m27, l28, m28, l29, m29, l30, m30,\
    l31, m31, l32, m32, l33, m33, l34, m34, l35, m35, l36, m36, l37, m37, l38, m38, l39, m39, l40,\
    m40, l41, m41, l42, m42, l43, m43, l44, m44, l45, m45, l46, m46, l47, m47, l48, m48, l49, m49,\
    l50, m50, l51, m51, l52, m52, l53, m53, l54, m54, l55, m55, l56, m56, l57, m57, l58, m58, l59,\
    m59, l60, m60, l61, m61, l62, m62, l63, m63, l64, m64, l65, m65, l66, m66, l67, m67, l68, m68,\
    l69, m69, l70, m70, l71, m71, l72, m72, l73, m73, l74, m74, l75, m75, l76, m76, l77, m77, l78,\
    m78, l79, m79, l80, m80, l81, m81, l82, m82, l83, m83, l84, m84, l85, m85, l86, m86, l87, m87,\
    l88, m88, l89, m89, l90, m90, l91, m91, l92, m92, l93, m93, l94, m94, l95, m95, l96, m96, l97,\
    m97, l98, m98, l99, m99, m100,g,q0, u0, q1, u1, q2, u2, q3, u3, q4, u4, q5, u5, q6, u6, q7, u7,\
    q8, u8, q9, u9, q10, u10, q11, u11,\
    q12, u12, q13, u13, q14, u14, q15, u15, q16, u16, q17, u17, q18, u18, q19, u19, q20, u20, q21,\
    u21, q22, u22, q23, u23, q24, u24, q25, u25, q26, u26, q27, u27, q28, u28, q29, u29, q30, u30,\
    q31, u31, q32, u32, q33, u33, q34, u34, q35, u35, q36, u36, q37, u37, q38, u38, q39, u39, q40,\
    u40, q41, u41, q42, u42, q43, u43, q44, u44, q45, u45, q46, u46, q47, u47, q48, u48, q49, u49,\
    q50, u50, q51, u51, q52, u52, q53, u53, q54, u54, q55, u55, q56, u56, q57, u57, q58, u58, q59,\
    u59, q60, u60, q61, u61, q62, u62, q63, u63, q64, u64, q65, u65, q66, u66, q67, u67, q68, u68,\
    q69, u69, q70, u70, q71, u71, q72, u72, q73, u73, q74, u74, q75, u75, q76, u76, q77, u77, q78,\
    u78, q79, u79, q80, u80, q81, u81, q82, u82, q83, u83, q84, u84, q85, u85, q86, u86, q87, u87,\
    q88, u88, q89, u89, q90, u90, q91, u91, q92, u92, q93, u93, q94, u94, q95, u95, q96, u96, q97,\
    u97, q98, u98, q99, u99, q100, u100, F, T1):
    import time
    m = 0
    for i in range(100):
        temps = time.time()
        _ = func(l0, m0, l1, m1, l2, m2, l3, m3, l4, m4, l5, m5, l6, m6, l7, m7, l8, m8, l9, m9,\
    l10, m10, l11, m11,\
    l12, m12, l13, m13, l14, m14, l15, m15, l16, m16, l17, m17, l18, m18, l19, m19, l20, m20, l21,\
    m21, l22, m22, l23, m23, l24, m24, l25, m25, l26, m26, l27, m27, l28, m28, l29, m29, l30, m30,\
    l31, m31, l32, m32, l33, m33, l34, m34, l35, m35, l36, m36, l37, m37, l38, m38, l39, m39, l40,\
    m40, l41, m41, l42, m42, l43, m43, l44, m44, l45, m45, l46, m46, l47, m47, l48, m48, l49, m49,\
    l50, m50, l51, m51, l52, m52, l53, m53, l54, m54, l55, m55, l56, m56, l57, m57, l58, m58, l59,\
    m59, l60, m60, l61, m61, l62, m62, l63, m63, l64, m64, l65, m65, l66, m66, l67, m67, l68, m68,\
    l69, m69, l70, m70, l71, m71, l72, m72, l73, m73, l74, m74, l75, m75, l76, m76, l77, m77, l78,\
    m78, l79, m79, l80, m80, l81, m81, l82, m82, l83, m83, l84, m84, l85, m85, l86, m86, l87, m87,\
    l88, m88, l89, m89, l90, m90, l91, m91, l92, m92, l93, m93, l94, m94, l95, m95, l96, m96, l97,\
    m97, l98, m98, l99, m99, m100,g,q0, u0, q1, u1, q2, u2, q3, u3, q4, u4, q5, u5, q6, u6, q7, u7,\
    q8, u8, q9, u9, q10, u10, q11, u11,\
    q12, u12, q13, u13, q14, u14, q15, u15, q16, u16, q17, u17, q18, u18, q19, u19, q20, u20, q21,\
    u21, q22, u22, q23, u23, q24, u24, q25, u25, q26, u26, q27, u27, q28, u28, q29, u29, q30, u30,\
    q31, u31, q32, u32, q33, u33, q34, u34, q35, u35, q36, u36, q37, u37, q38, u38, q39, u39, q40,\
    u40, q41, u41, q42, u42, q43, u43, q44, u44, q45, u45, q46, u46, q47, u47, q48, u48, q49, u49,\
    q50, u50, q51, u51, q52, u52, q53, u53, q54, u54, q55, u55, q56, u56, q57, u57, q58, u58, q59,\
    u59, q60, u60, q61, u61, q62, u62, q63, u63, q64, u64, q65, u65, q66, u66, q67, u67, q68, u68,\
    q69, u69, q70, u70, q71, u71, q72, u72, q73, u73, q74, u74, q75, u75, q76, u76, q77, u77, q78,\
    u78, q79, u79, q80, u80, q81, u81, q82, u82, q83, u83, q84, u84, q85, u85, q86, u86, q87, u87,\
    q88, u88, q89, u89, q90, u90, q91, u91, q92, u92, q93, u93, q94, u94, q95, u95, q96, u96, q97,\
    u97, q98, u98, q99, u99, q100, u100, F, T1)
        m += time.time()-temps
    print(m/100)


l0, m0, l1, m1, l2, m2, l3, m3, l4, m4, l5, m5, l6, m6, l7, m7, l8, m8, l9, m9, l10, m10, l11, m11,\
    l12, m12, l13, m13, l14, m14, l15, m15, l16, m16, l17, m17, l18, m18, l19, m19, l20, m20, l21,\
    m21, l22, m22, l23, m23, l24, m24, l25, m25, l26, m26, l27, m27, l28, m28, l29, m29, l30, m30,\
    l31, m31, l32, m32, l33, m33, l34, m34, l35, m35, l36, m36, l37, m37, l38, m38, l39, m39, l40,\
    m40, l41, m41, l42, m42, l43, m43, l44, m44, l45, m45, l46, m46, l47, m47, l48, m48, l49, m49,\
    l50, m50, l51, m51, l52, m52, l53, m53, l54, m54, l55, m55, l56, m56, l57, m57, l58, m58, l59,\
    m59, l60, m60, l61, m61, l62, m62, l63, m63, l64, m64, l65, m65, l66, m66, l67, m67, l68, m68,\
    l69, m69, l70, m70, l71, m71, l72, m72, l73, m73, l74, m74, l75, m75, l76, m76, l77, m77, l78,\
    m78, l79, m79, l80, m80, l81, m81, l82, m82, l83, m83, l84, m84, l85, m85, l86, m86, l87, m87,\
    l88, m88, l89, m89, l90, m90, l91, m91, l92, m92, l93, m93, l94, m94, l95, m95, l96, m96, l97,\
    m97, l98, m98, l99, m99, m100 =\
    sm.symbols("l0 m0 l1 m1 l2 m2 l3 m3 l4 m4 l5 m5 l6 m6 l7 m7 l8 m8 l9 m9 l10 m10 l11 m11 l12 m12 l13 m13 l14 m14 l15 m15 l16 m16 l17 m17 l18 m18 l19 m19 l20 m20 l21 m21 l22 m22 l23 m23 l24 m24 l25 m25 l26 m26 l27 m27 l28 m28 l29 m29 l30 m30 l31 m31 l32 m32 l33 m33 l34 m34 l35 m35 l36 m36 l37 m37 l38 m38 l39 m39 l40 m40 l41 m41 l42 m42 l43 m43 l44 m44 l45 m45 l46 m46 l47 m47 l48 m48 l49 m49 l50 m50 l51 m51 l52 m52 l53 m53 l54 m54 l55 m55 l56 m56 l57 m57 l58 m58 l59 m59 l60 m60 l61 m61 l62 m62 l63 m63 l64 m64 l65 m65 l66 m66 l67 m67 l68 m68 l69 m69 l70 m70 l71 m71 l72 m72 l73 m73 l74 m74 l75 m75 l76 m76 l77 m77 l78 m78 l79 m79 l80 m80 l81 m81 l82 m82 l83 m83 l84 m84 l85 m85 l86 m86 l87 m87 l88 m88 l89 m89 l90 m90 l91 m91 l92 m92 l93 m93 l94 m94 l95 m95 l96 m96 l97 m97 l98 m98 l99 m99 m100")

g = sm.symbols("g")

q0, q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, q14, q15, q16, q17, q18, q19, \
    q20, q21, q22, q23, q24, q25, q26, q27, q28, q29, q30, q31, q32, q33, q34, q35, q36,\
    q37, q38, q39, q40, q41, q42, q43, q44, q45, q46, q47, q48, q49, q50, q51, q52, q53, q54,\
    q55, q56, q57, q58, q59, q60, q61, q62, q63, q64, q65, q66, q67, q68, q69, q70, q71, q72,\
    q73, q74, q75, q76, q77, q78, q79, q80, q81, q82, q83, q84, q85, q86, q87, q88, q89, q90,\
    q91, q92, q93, q94, q95, q96, q97, q98, q99, q100=\
sm.Function('q0')(m0), sm.Function('q1')(m1), sm.Function('q2')(m2), sm.Function('q3')(m3), \
                       sm.Function('q4')(m4), sm.Function('q5')(m5), sm.Function('q6')(m6),\
                       sm.Function('q7')(m7), sm.Function('q8')(m8), sm.Function('q9')(m9),\
                       sm.Function('q10')(m10), sm.Function('q11')(m11), sm.Function('q12')(m12),\
                       sm.Function('q13')(m13), sm.Function('q14')(m14), sm.Function('q15')(m15),\
                       sm.Function('q16')(m16), sm.Function('q17')(m17), sm.Function('q18')(m18),\
                       sm.Function('q19')(m19), sm.Function('q20')(m20), sm.Function('q21')(m21),\
                       sm.Function('q22')(m22), sm.Function('q23')(m23), sm.Function('q24')(m24),\
                       sm.Function('q25')(m25), sm.Function('q26')(m26), sm.Function('q27')(m27),\
                       sm.Function('q28')(m28), sm.Function('q29')(m29), sm.Function('q30')(m30),\
                       sm.Function('q31')(m31), sm.Function('q32')(m32), sm.Function('q33')(m33),\
                       sm.Function('q34')(m34), sm.Function('q35')(m35), sm.Function('q36')(m36),\
                       sm.Function('q37')(m37), sm.Function('q38')(m38), sm.Function('q39')(m39),\
                       sm.Function('q40')(m40), sm.Function('q41')(m41), sm.Function('q42')(m42),\
                       sm.Function('q43')(m43), sm.Function('q44')(m44), sm.Function('q45')(m45),\
                       sm.Function('q46')(m46), sm.Function('q47')(m47), sm.Function('q48')(m48),\
                       sm.Function('q49')(m49), sm.Function('q50')(m50), sm.Function('q51')(m51),\
                       sm.Function('q52')(m52), sm.Function('q53')(m53), sm.Function('q54')(m54),\
                       sm.Function('q55')(m55), sm.Function('q56')(m56), sm.Function('q57')(m57),\
                       sm.Function('q58')(m58), sm.Function('q59')(m59), sm.Function('q60')(m60),\
                       sm.Function('q61')(m61), sm.Function('q62')(m62), sm.Function('q63')(m63),\
                       sm.Function('q64')(m64), sm.Function('q65')(m65), sm.Function('q66')(m66),\
                       sm.Function('q67')(m67), sm.Function('q68')(m68), sm.Function('q69')(m69), \
                       sm.Function('q70')(m70), sm.Function('q71')(m71), sm.Function('q72')(m72),\
                       sm.Function('q73')(m73), sm.Function('q74')(m74), sm.Function('q75')(m75),\
                       sm.Function('q76')(m76), sm.Function('q77')(m77), sm.Function('q78')(m78),\
                       sm.Function('q79')(m79), sm.Function('q80')(m80), sm.Function('q81')(m81),\
                       sm.Function('q82')(m82), sm.Function('q83')(m83), sm.Function('q84')(m84),\
                       sm.Function('q85')(m85), sm.Function('q86')(m86), sm.Function('q87')(m87),\
                       sm.Function('q88')(m88), sm.Function('q89')(m89), sm.Function('q90')(m90),\
                       sm.Function('q91')(m91), sm.Function('q92')(m92), sm.Function('q93')(m93),\
                       sm.Function('q94')(m94), sm.Function('q95')(m95), sm.Function('q96')(m96),\
                       sm.Function('q97')(m97), sm.Function('q98')(m98), sm.Function('q99')(m99), \
                       sm.Function('q100')(m100)

u0, u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15, u16, u17, u18, u19,\
    u20, u21, u22, u23, u24, u25, u26, u27, u28, u29, u30, u31, u32, u33, u34, u35, u36,\
    u37, u38, u39, u40, u41, u42, u43, u44, u45, u46, u47, u48, u49, u50, u51, u52, u53,\
    u54, u55, u56, u57, u58, u59, u60, u61, u62, u63, u64, u65, u66, u67, u68, u69, u70,\
    u71, u72, u73, u74, u75, u76, u77, u78, u79, u80, u81, u82, u83, u84, u85, u86, u87,\
    u88, u89, u90, u91, u92, u93, u94, u95, u96, u97, u98, u99, u100 =\
q0.diff(m0), q1.diff(m1), q2.diff(m2), q3.diff(m3), q4.diff(m4), q5.diff(m5), q6.diff(m6),\
q7.diff(m7), q8.diff(m8), q9.diff(m9), q10.diff(m10), q11.diff(m11), q12.diff(m12), q13.diff(m13),\
q14.diff(m14), q15.diff(m15), q16.diff(m16), q17.diff(m17), q18.diff(m18), q19.diff(m19),\
q20.diff(m20), q21.diff(m21), q22.diff(m22), q23.diff(m23), q24.diff(m24), q25.diff(m25),\
q26.diff(m26), q27.diff(m27), q28.diff(m28), q29.diff(m29), q30.diff(m30), q31.diff(m31),\
q32.diff(m32), q33.diff(m33), q34.diff(m34), q35.diff(m35), q36.diff(m36), q37.diff(m37),\
q38.diff(m38), q39.diff(m39), q40.diff(m40), q41.diff(m41), q42.diff(m42), q43.diff(m43), \
q44.diff(m44), q45.diff(m45), q46.diff(m46), q47.diff(m47), q48.diff(m48), q49.diff(m49),\
q50.diff(m50), q51.diff(m51), q52.diff(m52), q53.diff(m53), q54.diff(m54), q55.diff(m55),\
q56.diff(m56), q57.diff(m57), q58.diff(m58), q59.diff(m59), q60.diff(m60), q61.diff(m61),\
q62.diff(m62), q63.diff(m63), q64.diff(m64), q65.diff(m65), q66.diff(m66), q67.diff(m67),\
q68.diff(m68), q69.diff(m69), q70.diff(m70), q71.diff(m71), q72.diff(m72), q73.diff(m73),\
q74.diff(m74), q75.diff(m75), q76.diff(m76), q77.diff(m77), q78.diff(m78), q79.diff(m79),\
q80.diff(m80), q81.diff(m81), q82.diff(m82), q83.diff(m83), q84.diff(m84), q85.diff(m85),\
q86.diff(m86), q87.diff(m87), q88.diff(m88), q89.diff(m89), q90.diff(m90), q91.diff(m91),\
q92.diff(m92), q93.diff(m93), q94.diff(m94), q95.diff(m95), q96.diff(m96), q97.diff(m97),\
q98.diff(m98), q99.diff(m99), q100.diff(m100)


F, T1 = dynamicsymbols("F T1")

#generate matrix

kane1 = kane_with_derivatives(100)

avg_timing(sm.lambdify((l0, m0, l1, m1, l2, m2, l3, m3, l4, m4, l5, m5, l6, m6, l7, m7, l8, m8, l9,\
                        m9, l10, m10, l11, m11,\
    l12, m12, l13, m13, l14, m14, l15, m15, l16, m16, l17, m17, l18, m18, l19, m19, l20, m20, l21,\
    m21, l22, m22, l23, m23, l24, m24, l25, m25, l26, m26, l27, m27, l28, m28, l29, m29, l30, m30,\
    l31, m31, l32, m32, l33, m33, l34, m34, l35, m35, l36, m36, l37, m37, l38, m38, l39, m39, l40,\
    m40, l41, m41, l42, m42, l43, m43, l44, m44, l45, m45, l46, m46, l47, m47, l48, m48, l49, m49,\
    l50, m50, l51, m51, l52, m52, l53, m53, l54, m54, l55, m55, l56, m56, l57, m57, l58, m58, l59,\
    m59, l60, m60, l61, m61, l62, m62, l63, m63, l64, m64, l65, m65, l66, m66, l67, m67, l68, m68,\
    l69, m69, l70, m70, l71, m71, l72, m72, l73, m73, l74, m74, l75, m75, l76, m76, l77, m77, l78,\
    m78, l79, m79, l80, m80, l81, m81, l82, m82, l83, m83, l84, m84, l85, m85, l86, m86, l87, m87,\
    l88, m88, l89, m89, l90, m90, l91, m91, l92, m92, l93, m93, l94, m94, l95, m95, l96, m96, l97,\
    m97, l98, m98, l99, m99, m100,g,q0, u0, q1, u1, q2, u2, q3, u3, q4, u4, q5, u5, q6, u6, q7, u7,\
                        q8, u8, q9, u9, q10, u10, q11, u11,\
    q12, u12, q13, u13, q14, u14, q15, u15, q16, u16, q17, u17, q18, u18, q19, u19, q20, u20, q21,\
    u21, q22, u22, q23, u23, q24, u24, q25, u25, q26, u26, q27, u27, q28, u28, q29, u29, q30, u30,\
    q31, u31, q32, u32, q33, u33, q34, u34, q35, u35, q36, u36, q37, u37, q38, u38, q39, u39, q40,\
    u40, q41, u41, q42, u42, q43, u43, q44, u44, q45, u45, q46, u46, q47, u47, q48, u48, q49, u49,\
    q50, u50, q51, u51, q52, u52, q53, u53, q54, u54, q55, u55, q56, u56, q57, u57, q58, u58, q59,\
    u59, q60, u60, q61, u61, q62, u62, q63, u63, q64, u64, q65, u65, q66, u66, q67, u67, q68, u68,\
    q69, u69, q70, u70, q71, u71, q72, u72, q73, u73, q74, u74, q75, u75, q76, u76, q77, u77, q78,\
    u78, q79, u79, q80, u80, q81, u81, q82, u82, q83, u83, q84, u84, q85, u85, q86, u86, q87, u87,\
    q88, u88, q89, u89, q90, u90, q91, u91, q92, u92, q93, u93, q94, u94, q95, u95, q96, u96, q97,\
    u97, q98, u98, q99, u99, q100, u100, F, T1 ), kane1),\
       40, 8, 19, 40, 4, 39, 29, 42, 22, 0, 40, 21, 6, 31, 30, 12, 25, 36, 15, 43, 48, 31, 8, 45, 45,\
           11, 5, 21, 0, 27, 14, 4, 36, 23, 17, 2, 22, 40, 26, 47, 21, 39, 32, 34, 1, 0, 39, 10, 33, 4,\
           37, 22, 47, 29, 41, 38, 0, 15, 20, 2, 47, 24, 15, 25, 14, 4, 37, 30, 17, 37, 2, 0, 25, 5, 14,\
           10, 26, 12, 12, 21, 49, 7, 30, 23, 40, 23, 48, 6, 24, 5, 43, 33, 22, 11, 43, 47, 0, 25, 15, 20,\
           12, 6, 41, 10, 21, 34, 36, 48, 4, 5, 44, 47, 10, 36, 37, 49, 49, 28, 23, 28, 40, 21, 42, 40, 45,\
           29, 14, 20, 3, 4, 14, 6, 20, 40, 48, 22, 34, 5, 5, 28, 47, 17, 1, 24, 21, 10, 26, 24, 28, 5, 34,\
           22, 30, 15, 20, 46, 11, 44, 36, 43, 14, 5, 10, 5, 38, 37, 26, 26, 28, 38, 20, 25, 42, 34, 14, 47,\
           2, 18, 9, 23, 6, 24, 16, 33, 6, 7, 44, 39, 40, 49, 47, 27, 27, 10, 17, 5, 44, 1, 10, 30, 39, 18,\
           27, 31, 42, 23, 8, 14, 48, 12, 46, 27, 39, 47, 39, 38, 19, 43, 48, 6, 49, 6, 37, 17, 34, 35, 19,\
           18, 9, 22, 19, 7, 48, 14, 42, 8, 48, 21, 9, 34, 17, 49, 17, 24, 36, 48, 34, 9, 20, 2, 30, 38, 31,\
           30, 29, 7, 19, 45, 49, 34, 44, 9, 43, 17, 6, 34, 34, 43, 45, 7, 8, 22, 11, 22, 24, 29, 7, 32, 13,\
           26, 36, 27, 35, 37, 25, 13, 48, 23, 47, 16, 14, 24, 42, 10, 21, 19, 27, 12, 37, 7, 22, 16, 1, 26,\
           46, 32, 26, 42, 19, 30, 6, 12, 49, 14, 45, 37, 19, 24, 23, 46, 47, 44, 12, 44, 20, 19, 4, 36, 20,\
           0, 16, 39, 4, 24, 24, 20, 43, 21, 14, 26, 34, 3, 48, 35, 46, 43, 44, 11, 18, 24, 6, 2, 48, 13, 8,\
           29, 31, 44, 15, 44, 21, 49, 47, 4, 17, 15, 46, 25, 1, 31, 5, 39, 44, 0, 31, 28, 2, 22, 17, 40, 47,\
           36, 17, 38, 42, 17, 21, 21, 26, 35, 26, 32, 45, 28, 26, 35, 14, 37, 24, 25, 13, 40, 43, 21, 35, 15)
#and for cse enabled

avg_timing(sm.lambdify((l0, m0, l1, m1, l2, m2, l3, m3, l4, m4, l5, m5, l6, m6, l7, m7, l8, m8, l9,\
                        m9, l10, m10, l11, m11,\
    l12, m12, l13, m13, l14, m14, l15, m15, l16, m16, l17, m17, l18, m18, l19, m19, l20, m20, l21,\
    m21, l22, m22, l23, m23, l24, m24, l25, m25, l26, m26, l27, m27, l28, m28, l29, m29, l30, m30,\
    l31, m31, l32, m32, l33, m33, l34, m34, l35, m35, l36, m36, l37, m37, l38, m38, l39, m39, l40,\
    m40, l41, m41, l42, m42, l43, m43, l44, m44, l45, m45, l46, m46, l47, m47, l48, m48, l49, m49,\
    l50, m50, l51, m51, l52, m52, l53, m53, l54, m54, l55, m55, l56, m56, l57, m57, l58, m58, l59,\
    m59, l60, m60, l61, m61, l62, m62, l63, m63, l64, m64, l65, m65, l66, m66, l67, m67, l68, m68,\
    l69, m69, l70, m70, l71, m71, l72, m72, l73, m73, l74, m74, l75, m75, l76, m76, l77, m77, l78,\
    m78, l79, m79, l80, m80, l81, m81, l82, m82, l83, m83, l84, m84, l85, m85, l86, m86, l87, m87,\
    l88, m88, l89, m89, l90, m90, l91, m91, l92, m92, l93, m93, l94, m94, l95, m95, l96, m96, l97,\
    m97, l98, m98, l99, m99, m100,g,q0, u0, q1, u1, q2, u2, q3, u3, q4, u4, q5, u5, q6, u6, q7, u7,\
                        q8, u8, q9, u9, q10, u10, q11, u11,\
    q12, u12, q13, u13, q14, u14, q15, u15, q16, u16, q17, u17, q18, u18, q19, u19, q20, u20, q21,\
    u21, q22, u22, q23, u23, q24, u24, q25, u25, q26, u26, q27, u27, q28, u28, q29, u29, q30, u30,\
    q31, u31, q32, u32, q33, u33, q34, u34, q35, u35, q36, u36, q37, u37, q38, u38, q39, u39, q40,\
    u40, q41, u41, q42, u42, q43, u43, q44, u44, q45, u45, q46, u46, q47, u47, q48, u48, q49, u49,\
    q50, u50, q51, u51, q52, u52, q53, u53, q54, u54, q55, u55, q56, u56, q57, u57, q58, u58, q59,\
    u59, q60, u60, q61, u61, q62, u62, q63, u63, q64, u64, q65, u65, q66, u66, q67, u67, q68, u68,\
    q69, u69, q70, u70, q71, u71, q72, u72, q73, u73, q74, u74, q75, u75, q76, u76, q77, u77, q78,\
    u78, q79, u79, q80, u80, q81, u81, q82, u82, q83, u83, q84, u84, q85, u85, q86, u86, q87, u87,\
    q88, u88, q89, u89, q90, u90, q91, u91, q92, u92, q93, u93, q94, u94, q95, u95, q96, u96, q97,\
    u97, q98, u98, q99, u99, q100, u100, F, T1 ), kane1, cse=True),\
       40, 8, 19, 40, 4, 39, 29, 42, 22, 0, 40, 21, 6, 31, 30, 12, 25, 36, 15, 43, 48, 31, 8, 45, 45,\
           11, 5, 21, 0, 27, 14, 4, 36, 23, 17, 2, 22, 40, 26, 47, 21, 39, 32, 34, 1, 0, 39, 10, 33, 4,\
           37, 22, 47, 29, 41, 38, 0, 15, 20, 2, 47, 24, 15, 25, 14, 4, 37, 30, 17, 37, 2, 0, 25, 5, 14,\
           10, 26, 12, 12, 21, 49, 7, 30, 23, 40, 23, 48, 6, 24, 5, 43, 33, 22, 11, 43, 47, 0, 25, 15, 20,\
           12, 6, 41, 10, 21, 34, 36, 48, 4, 5, 44, 47, 10, 36, 37, 49, 49, 28, 23, 28, 40, 21, 42, 40, 45,\
           29, 14, 20, 3, 4, 14, 6, 20, 40, 48, 22, 34, 5, 5, 28, 47, 17, 1, 24, 21, 10, 26, 24, 28, 5, 34,\
           22, 30, 15, 20, 46, 11, 44, 36, 43, 14, 5, 10, 5, 38, 37, 26, 26, 28, 38, 20, 25, 42, 34, 14, 47,\
           2, 18, 9, 23, 6, 24, 16, 33, 6, 7, 44, 39, 40, 49, 47, 27, 27, 10, 17, 5, 44, 1, 10, 30, 39, 18,\
           27, 31, 42, 23, 8, 14, 48, 12, 46, 27, 39, 47, 39, 38, 19, 43, 48, 6, 49, 6, 37, 17, 34, 35, 19,\
           18, 9, 22, 19, 7, 48, 14, 42, 8, 48, 21, 9, 34, 17, 49, 17, 24, 36, 48, 34, 9, 20, 2, 30, 38, 31,\
           30, 29, 7, 19, 45, 49, 34, 44, 9, 43, 17, 6, 34, 34, 43, 45, 7, 8, 22, 11, 22, 24, 29, 7, 32, 13,\
           26, 36, 27, 35, 37, 25, 13, 48, 23, 47, 16, 14, 24, 42, 10, 21, 19, 27, 12, 37, 7, 22, 16, 1, 26,\
           46, 32, 26, 42, 19, 30, 6, 12, 49, 14, 45, 37, 19, 24, 23, 46, 47, 44, 12, 44, 20, 19, 4, 36, 20,\
           0, 16, 39, 4, 24, 24, 20, 43, 21, 14, 26, 34, 3, 48, 35, 46, 43, 44, 11, 18, 24, 6, 2, 48, 13, 8,\
           29, 31, 44, 15, 44, 21, 49, 47, 4, 17, 15, 46, 25, 1, 31, 5, 39, 44, 0, 31, 28, 2, 22, 17, 40, 47,\
           36, 17, 38, 42, 17, 21, 21, 26, 35, 26, 32, 45, 28, 26, 35, 14, 37, 24, 25, 13, 40, 43, 21, 35, 15)
