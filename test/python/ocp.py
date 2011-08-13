from casadi import *
import casadi as c
from numpy import *
import unittest
from types import *
from helpers import *

class OCPtests(casadiTestCase):
  def testdiscrete(self):
    self.message("Linear-quadratic problem, discrete, using IPOPT")
    # inspired by www.cs.umsl.edu/~janikow/publications/1992/GAforOpt/text.pdf
    a=1.0
    b=1.0
    q=1.0
    s=1.0
    r=1.0
    x0=100
    
    N=100
    
    X=symbolic("X",N+1)
    U=symbolic("U",N)
    
    V = vertcat([X,U])
    
    cost = 0
    for i in range(N):
      cost = cost + s*X[i]**2+r*U[i]**2
    cost = cost + q*X[N]**2
    
    f = SXFunction([V],[cost])
    
    g = SXFunction([V],[vertcat([X[0]-x0,X[1:,0]-(a*X[:N,0]+b*U)])])
    
    solver = IpoptSolver(f,g)
    solver.setOption("tol",1e-5)
    solver.setOption("hessian_approximation", "limited-memory")
    solver.setOption("max_iter",100)
    solver.setOption("print_level",0)
    solver.init()
    solver.input(NLP_LBX).set([-1000 for i in range(V.size())])
    solver.input(NLP_UBX).set([1000 for i in range(V.size())])
    solver.input(NLP_LBG).set([0 for i in range(N+1)])
    solver.input(NLP_UBG).set([0 for i in range(N+1)])
    solver.solve()
    ocp_sol=solver.output(NLP_COST)[0]
    # solve the ricatti equation exactly
    K = q+0.0
    for i in range(N):
      K = s+r*a**2*K/(r+b**2*K)
    exact_sol=K * x0**2
    self.assertAlmostEqual(ocp_sol,exact_sol,10,"Linear-quadratic problem solution using IPOPT")
    
  def test_singleshooting(self):
    self.message("Single shooting")
    p0 = 0.2
    y0= 1
    yc0=dy0=0
    te=0.4

    t=symbolic("t")
    q=symbolic("y",2,1)
    p=symbolic("p",1,1)
    # y
    # y'
    f=SXFunction([t,q,p,[]],[vertcat([q[1],p[0]+q[1]**2 ])])
    f.init()
    
    integrator = CVodesIntegrator(f)
    integrator.setOption("reltol",1e-15)
    integrator.setOption("abstol",1e-15)
    integrator.setOption("verbose",True)
    integrator.setOption("steps_per_checkpoint",10000)
    integrator.setOption("t0",0)
    integrator.setOption("tf",te)

    integrator.init()

    var = MX("var",2,1)
    par = MX("par",1,1)
    parMX= par
    
    q0   = vertcat([var[0],par])
    par  = var[1]
    qend=integrator([q0,par,MX(2,1)])
    
    parc = MX(0)
    
    f = MXFunction([var,parMX],[qend[0]])
    f.init()
    fc = MXFunction([var],[-f([var,parc])])
    fc.init()
    solver = IpoptSolver(fc)
    solver.setOption("tol",1e-12)
    solver.setOption("hessian_approximation", "limited-memory")
    solver.setOption("max_iter",10)
    solver.setOption("derivative_test","first-order")
    solver.setOption("print_level",0)
    solver.init()
    solver.input(NLP_LBX).set([-1, -1])
    solver.input(NLP_UBX).set([1, 0.2])
    solver.solve()
    self.assertAlmostEqual(solver.output(NLP_X_OPT)[0],1,8,"X_opt")
    self.assertAlmostEqual(solver.output(NLP_X_OPT)[1],0.2,8,"X_opt")
    
    self.assertAlmostEqual(solver.output(NLP_LAMBDA_UBX)[0],1,8,"Cost should be linear in y0")
    self.assertAlmostEqual(solver.output(NLP_LAMBDA_UBX)[1],(sqrt(p0)*(te*yc0**2-yc0+p0*te)*tan(arctan(yc0/sqrt(p0))+sqrt(p0)*te)+yc0**2)/(2*p0*yc0**2+2*p0**2),8,"Cost should be linear in y0")
    self.assertAlmostEqual(-solver.output(NLP_COST)[0],(2*y0-log(yc0**2/p0+1))/2-log(cos(arctan(yc0/sqrt(p0))+sqrt(p0)*te)),7,"Cost")
    self.assertAlmostEqual(solver.output(NLP_LAMBDA_LBX)[0],0,8,"Constraint is supposed to be unactive")
    self.assertAlmostEqual(solver.output(NLP_LAMBDA_LBX)[1],0,8,"Constraint is supposed to be unactive")
  
  def test_singleshooting2(self):
    self.message("Single shooting 2")
    p0 = 0.2
    y0= 0.2
    yc0=dy0=0.1
    te=0.4

    t=symbolic("t")
    q=symbolic("y",2,1)
    p=symbolic("p",1,1)
    # y
    # y'
    f=SXFunction([t,q,p,[]],[vertcat([q[1],p[0]+q[1]**2 ])])
    f.init()
    
    integrator = CVodesIntegrator(f)
    integrator.setOption("reltol",1e-15)
    integrator.setOption("abstol",1e-15)
    integrator.setOption("verbose",True)
    integrator.setOption("steps_per_checkpoint",10000)
    integrator.setOption("t0",0)
    integrator.setOption("tf",te)

    integrator.init()

    var = MX("var",2,1)
    par = MX("par",1,1)
    
    q0   = vertcat([var[0],par])
    parl  = var[1]
    qend=integrator([q0,parl,MX(2,1)])
    
    parc = MX(dy0)
    
    f = MXFunction([var,par],[qend[0]])
    f.init()
    fc = MXFunction([var],[-f([var,parc])])
    fc.init()
    
    g = MXFunction([var],[var[0]-var[1]])
    g.init()
    
    solver = IpoptSolver(fc,g)
    solver.setOption("tol",1e-12)
    solver.setOption("hessian_approximation", "limited-memory")
    solver.setOption("max_iter",10)
    solver.setOption("derivative_test","first-order")
    #solver.setOption("print_level",0)
    solver.init()
    solver.input(NLP_LBX).set([-1, -1])
    solver.input(NLP_UBX).set([1, 0.2])
    solver.input(NLP_LBG).set([-1])
    solver.input(NLP_UBG).set([0])
    solver.solve()

    self.assertAlmostEqual(solver.output(NLP_X_OPT)[0],0.2,6,"X_opt")
    self.assertAlmostEqual(solver.output(NLP_X_OPT)[1],0.2,6,"X_opt")
    
    self.assertAlmostEqual(solver.output(NLP_LAMBDA_UBX)[0],0,8,"Constraint is supposed to be unactive")
    dfdp0 = (sqrt(p0)*(te*yc0**2-yc0+p0*te)*tan(arctan(yc0/sqrt(p0))+sqrt(p0)*te)+yc0**2)/(2*p0*yc0**2+2*p0**2)
    self.assertAlmostEqual(solver.output(NLP_LAMBDA_UBX)[1],1+dfdp0,8)
    self.assertAlmostEqual(solver.output(NLP_LAMBDA_OPT)[0],1,8)
    self.assertAlmostEqual(-solver.output(NLP_COST)[0],(2*y0-log(yc0**2/p0+1))/2-log(cos(arctan(yc0/sqrt(p0))+sqrt(p0)*te)),7,"Cost")
    self.assertAlmostEqual(solver.output(NLP_LAMBDA_LBX)[0],0,8,"Constraint is supposed to be unactive")
    self.assertAlmostEqual(solver.output(NLP_LAMBDA_LBX)[1],0,8,"Constraint is supposed to be unactive") 
    
  def test_XML(self):
    self.message("JModelica XML parsing")
    parser = FMIParser('data/cstr.xml')
    ocp = parser.parse()
    
    self.assertEqual(ocp.t0,0)
    self.assertEqual(ocp.tf,150)
    #self.assertFalse(ocp.t0_free)
    #self.assertFalse(ocp.tf_free)
    self.assertTrue(len(ocp.lterm)==0)
    self.assertTrue(len(ocp.mterm)==1)
    m = ocp.mterm[0]
    self.assertTrue(isinstance(m,SX))
    self.assertTrue(isinstance(ocp.t_,SX))
    self.assertEquals(str(m),'cost.atTime(150)')
    print dir(ocp)
    self.assertEquals(len(ocp.implicit_fcn_),3)
    self.assertEquals(len(ocp.implicit_fcn_),3)
    self.assertEquals(len(ocp.x_),3) # there are three states
    (c,T,cost) = tuple(ocp.x_)
    self.assertTrue(isinstance(c,Variable))
    self.assertEquals(c.getName(),"cstr.c")
    self.assertEquals(T.getName(),"cstr.T")
    self.assertEquals(cost.getName(),"cost")
    self.assertEquals(c.getNominal(),1000)
    
    w=ocp.variables_.subByName('cstr')
    #self.assertEquals(w("T"),T)
    self.assertEquals(w.subByName("T").var_.getName(),"cstr.T")
   
    #print c.atTime(0)
       
    c = c.var()
    T = T.var()
    cost = cost.var()
    
    self.assertTrue(w.subByName("T").var_.var().isEqual(T)) 
        
    u = ocp.u_[0].var()
    self.assertEquals(len(ocp.path_fcn_),3)
    #self.assertEquals(len(ocp.cfcn_lb),3)
    #self.assertEquals(len(ocp.cfcn_ub),3)
    #self.assertTrue(ocp.cfcn[0].isEqual(T)) 
    #self.assertTrue(ocp.cfcn[1].isEqual(u)) 
    #self.assertTrue(ocp.cfcn[2].isEqual(u)) 
    #self.assertTrue(ocp.cfcn_lb[0].isMinusInf()) 
    #self.assertEquals(ocp.cfcn_lb[1].getValue(),230) 
    #self.assertTrue(ocp.cfcn_lb[2].isMinusInf()) 
    #self.assertEquals(ocp.cfcn_ub[0].getValue(),350) 
    #self.assertTrue(ocp.cfcn_ub[1].isInf())
    #self.assertEquals(ocp.cfcn_ub[2].getValue(),370) 
    print ocp.initial_eq_
    print c,T,cost
    #print c.atTime(0)
    f=SXFunction([[c,T,cost]],[ocp.initial_eq_])
    f.init()
    return 
    f.evaluate()
    self.checkarray(f.output(),matrix([-956.271065,-250.051971,0]).T,"initeq")

    
    mystates = []

  def testMSclass_prim(self):
    self.message("CasADi multiple shooting class")
    
    ns = 20
    nx = 3
    nu = 2
    np = 0
    nh = 0
    tf = 0.2
    
    x0 = symbolic("x0",nx)
    p = symbolic("p",nu)
    xp0 = symbolic("x0",nx)
    xf = x0 + p[0]
    dynamics = SXFunction({'NUM': INTEGRATOR_NUM_IN, INTEGRATOR_X0: x0,INTEGRATOR_P: p, INTEGRATOR_XP0: xp0},{'NUM': INTEGRATOR_NUM_OUT, INTEGRATOR_XF: xf,INTEGRATOR_XPF: xp0})
    mayer = SXFunction([x0],[7*x0[0]])
    ms = MultipleShooting(dynamics,mayer)
    ms.setOption("number_of_grid_points",ns)
    ms.setOption("final_time",tf)
    ms.setOption("nlp_solver",IpoptSolver)
    ms.init()
    self.checkarray(linspace(0,tf,ns+1),ms.input(OCP_T),"timegrid")
    
    for i in [OCP_LBX,OCP_UBX,OCP_X_INIT]:
      self.checkarray(ms.input(i).shape,(nx,ns+1),"shape")
      

    self.checkarray(ms.input(OCP_LBXP).shape,(nx,ns+1),"shape")
    self.checkarray(ms.input(OCP_UBXP).shape,(nx,ns+1),"shape")
    self.checkarray(ms.input(OCP_XP_INIT).shape,(0,0),"shape") # is this a bug?
    
    for i in [OCP_LBU,OCP_UBU,OCP_U_INIT]:
      self.checkarray(ms.input(i).shape,(nu,ns),"shape")
    
    for i in [OCP_LBP,OCP_UBP,OCP_P_INIT]:
      self.checkarray(ms.input(i).shape,(np,1),"shape")

    for i in [OCP_LBH,OCP_UBH]:
      self.checkarray(ms.input(i).shape,(nh,ns+1),"shape")
      
    ns = 20
    nx = 3
    nu = 2
    np = 4
    nh = 2
    tf = 0.2
    
    x0 = symbolic("x0",nx)
    p = symbolic("p",nu+np)
    xp0 = symbolic("x0",nx)
    xf = x0 + p[0]
    dynamics = SXFunction({INTEGRATOR_X0: x0,INTEGRATOR_P: p, INTEGRATOR_XP0: xp0},{INTEGRATOR_XF: xf,INTEGRATOR_XPF: xp0})
    mayer = SXFunction([x0],[7*x0[0]])
    
    t = SX("t")
    cfcn = SXFunction({DAE_T : t, DAE_Y: x0, DAE_P: p, DAE_YDOT: xp0},[x0[:nh,0]])
    cfcn.init()
    
    ms = MultipleShooting(dynamics,mayer,cfcn)
    ms.setOption("number_of_grid_points",ns)
    ms.setOption("number_of_parameters",np)
    ms.setOption("final_time",tf)
    ms.setOption("nlp_solver",IpoptSolver)
    ms.init()
    self.checkarray(linspace(0,tf,ns+1),ms.input(OCP_T),"timegrid")
    
    for i in [OCP_LBX,OCP_UBX,OCP_X_INIT]:
      self.checkarray(ms.input(i).shape,(nx,ns+1),"shape")
      

    self.checkarray(ms.input(OCP_LBXP).shape,(nx,ns+1),"shape")
    self.checkarray(ms.input(OCP_UBXP).shape,(nx,ns+1),"shape")
    self.checkarray(ms.input(OCP_XP_INIT).shape,(0,0),"shape") # is this a bug?
    
    for i in [OCP_LBU,OCP_UBU,OCP_U_INIT]:
      self.checkarray(ms.input(i).shape,(nu,ns),"shape")
    
    for i in [OCP_LBP,OCP_UBP,OCP_P_INIT]:
      self.checkarray(ms.input(i).shape,(np,1),"shape")

    for i in [OCP_LBH,OCP_UBH]:
      self.checkarray(ms.input(i).shape,(nh,ns+1),"shape")

  def testMSclassSimple(self):
    self.message("CasADi multiple shooting class: simple example")
    """
    The problem consists of a harmonic oscilator and a power harvester.
    
    max     int_0^T u(t) * x(t)
     u,x,y
            s.t    x'(t) = y(t)
                   y'(t) = -x(t)
                   -1  <=   u(t) <= 1
                   x(0) = 1
                   y(0) = 0
                   
    The trivial solution for u(t) is u(t)=sign(x(t))
    
    
    """
    te = 2*pi
    N = 20
    t=SX("t")
    y=symbolic("y",3,1)
    yd=symbolic("yd",3,1)
    p=SX("p")
    f=SXFunction({'NUM':DAE_NUM_IN, DAE_T: t, DAE_Y: y, DAE_P: p,DAE_YDOT: yd},[[y[1,0],-y[0,0],p*y[0,0]]])
    f.init()
    
    integrator = CVodesIntegrator(f)
    integrator.setOption("reltol",1e-9)
    integrator.setOption("abstol",1e-9)
    integrator.setOption("steps_per_checkpoint",10000)
    integrator.setOption("t0",0)
    integrator.setOption("tf",te/N)
    integrator.init()
    
    mayer = SXFunction([y],[-y[2]])
    mayer.init()
    
    ms = MultipleShooting(integrator,mayer)
    ms.setOption("number_of_grid_points",N);
    ms.setOption("final_time",te);
    
    ms.setOption("nlp_solver",IpoptSolver)
    nlp_solver_options = {}
    nlp_solver_options["tol"] = 1e-10
    nlp_solver_options["hessian_approximation"] = "limited-memory"
    nlp_solver_options["max_iter"] = 100
    nlp_solver_options["linear_solver"] = "ma57"
    nlp_solver_options["derivative_test"] = "first-order"
    ms.setOption("nlp_solver_options",nlp_solver_options)
    ms.init()
    
    ms.input(OCP_LBX).setAll(-inf)
    ms.input(OCP_UBX).setAll(inf)
    ms.input(OCP_X_INIT).setAll(0)
    
    ms.input(OCP_LBU).setAll(-1)
    ms.input(OCP_UBU).setAll(1)
    ms.input(OCP_U_INIT).setAll(0)
    
    ms.input(OCP_LBX)[0,0] = 1
    ms.input(OCP_UBX)[0,0] = 1
   
    ms.input(OCP_LBX)[1,0] = 0
    ms.input(OCP_UBX)[1,0] = 0
 
    ms.input(OCP_LBX)[2,0] = 0
    ms.input(OCP_UBX)[2,0] = 0
    
    ms.solve()
    
    self.checkarray(sign(matrix(ms.output(OCP_X_OPT))[0,:-1]),ms.output(OCP_U_OPT),"solution")


  def testMSclassSimple2(self):
    return 
    self.message("CasADi multiple shooting class: simple example 2")
    """
    The problem consists of a harmonic oscilator and a power harvester.
    
    max     x(tend)
     u,x,a
            s.t    x'(t) = a*x(t) + u(t)
                   -1  <=   u(t) <= 1
                   -2  <=   a <= 2
                   x(0) = 0
                   
    The trivial solution for u(t) is u(t)=1, a = 1
    x(t) = 1/2 (e^(2 t)-1)
    
    
    """
    te = 1
    N = 20
    t=SX("t")
    x=SX("x")
    xd=SX("xd")
    a=SX("a")
    u=SX("u")
    f=SXFunction({'NUM':DAE_NUM_IN, DAE_T: t, DAE_Y: x, DAE_P: [a,u],DAE_YDOT: xd},[[a*x+u]])
    f.init()
    
    integrator = CVodesIntegrator(f)
    integrator.setOption("reltol",1e-9)
    integrator.setOption("abstol",1e-9)
    integrator.setOption("steps_per_checkpoint",10000)
    integrator.setOption("t0",0)
    integrator.setOption("tf",te/N)
    integrator.init()
    
    mayer = SXFunction([x],[-x])
    mayer.init()
    
    ms = MultipleShooting(integrator,mayer)
    ms.setOption("number_of_grid_points",N);
    ms.setOption("number_of_parameters",1);
    ms.setOption("final_time",te);

    ms.setOption("nlp_solver",IpoptSolver)
    nlp_solver_options = {}
    nlp_solver_options["tol"] = 1e-10
    nlp_solver_options["hessian_approximation"] = "limited-memory"
    nlp_solver_options["max_iter"] = 100
    nlp_solver_options["linear_solver"] = "ma57"
    nlp_solver_options["derivative_test"] = "first-order"
    ms.setOption("nlp_solver_options",nlp_solver_options)
    ms.init()
    
    ms.input(OCP_LBX).setAll(-inf)
    ms.input(OCP_UBX).setAll(inf)
    ms.input(OCP_X_INIT).setAll(0)
    
    ms.input(OCP_LBU).setAll(-1)
    ms.input(OCP_UBU).setAll(1)
    ms.input(OCP_U_INIT).setAll(0)
    
    ms.input(OCP_LBP).setAll(-2)
    ms.input(OCP_UBP).setAll(2)
    ms.input(OCP_P_INIT).setAll(0)
    
    ms.input(OCP_LBX)[0,0] = 0
    ms.input(OCP_UBX)[0,0] = 0
   
    ms.solve()
    
  
if __name__ == '__main__':
    unittest.main()

