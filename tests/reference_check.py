"""Independent SciPy check; not execution of the MATLAB source."""
import json
import numpy as np
from scipy.sparse import lil_matrix, coo_matrix
from scipy.sparse.linalg import spsolve

def solve(eps, x, y, left, right, bottom, top, source=None, initial=None, tol=1e-10):
    nx, ny = len(x), len(y)
    h, p = x[1]-x[0], y[1]-y[0]
    rx, ry = eps/h**2, 1/p**2
    if source is None: source = np.zeros((ny,nx))
    a = lil_matrix((nx*ny,nx*ny)); b=source.copy()
    active=[]
    for j in range(ny):
        for i in range(nx):
            q = j*nx+i
            if i in (0,nx-1):
                a[q,q]=1; b[j,i]=left[j] if i==0 else right[j]
            else:
                active.append(q)
                a[q,q]=2*rx+2*ry; a[q,q-1]=-rx; a[q,q+1]=-rx
                if j==0: a[q,q+nx]=-2*ry; b[j,i]-=2*bottom[i]/p
                elif j==ny-1: a[q,q-nx]=-2*ry; b[j,i]+=2*top[i]/p
                else: a[q,q-nx]=-ry; a[q,q+nx]=-ry
    a=a.tocsr(); b=b.ravel(); active=np.array(active)
    def convection(u):
        c=u[active]/(2*h)
        return coo_matrix((np.r_[-c,c],(np.r_[active,active],np.r_[active-1,active+1])),shape=a.shape).tocsr()
    def residual(u):
        # Independently evaluate the physical ghost-value stencil.
        z=u.reshape(ny,nx)
        ext=np.vstack((z[1]-2*p*bottom,z,z[-2]+2*p*top))
        r=z[:,1:-1]*(z[:,2:]-z[:,:-2])/(2*h)
        r-=(ext[2:,1:-1]-2*z[:,1:-1]+ext[:-2,1:-1])/p**2
        r-=eps*(z[:,2:]-2*z[:,1:-1]+z[:,:-2])/h**2
        r-=source[:,1:-1]
        return max(np.max(np.abs(r)),np.max(np.abs(z[:,0]-left)),np.max(np.abs(z[:,-1]-right)))
    u=spsolve(a,b) if initial is None else initial.ravel().copy()
    u[::nx]=left;u[nx-1::nx]=right
    res=np.max(np.abs((a+convection(u))@u-b))
    history=[]
    if res <= tol: return u.reshape(ny,nx),history,res,residual(u)
    for k in range(300):
        v=spsolve(a+convection(u),b)
        step=1
        while True:
            trial=u+step*(v-u)
            r=np.max(np.abs((a+convection(trial))@trial-b))
            if r <= res*(1-1e-4*step) or r <= tol: break
            if step <= 1/1024: raise RuntimeError(('stalled',k,res,r))
            step=max(step/2,1/1024)
        change=np.max(abs(trial-u))/max(1,np.max(abs(u)))
        u=trial;res=r;history.append([k+1,change,res,step])
        if change <= 1e-11 and res<=tol: return u.reshape(ny,nx),history,res,residual(u)
    raise RuntimeError(('limit',res))

x=np.linspace(0,1,13);y=np.linspace(0,1,11);X,Y=np.meshgrid(x,y)
z=np.zeros_like(X)
u,h,r,ind=solve(.4,x,y,z[:,0],z[:,-1],np.zeros(len(x)),np.zeros(len(x)))
assert np.max(abs(u))<1e-12
exact=Y**2+.3*Y+.2*X+.1
u,h,r,ind=solve(.4,x,y,exact[:,0],exact[:,-1],np.full(len(x),.3),np.full(len(x),2.3),.2*exact-2)
quad={'max_error':float(np.max(abs(u-exact))),'iterations':len(h),'matrix_residual':r,'stencil_residual':ind}
assert quad['max_error']<1e-8
refine=[]
for n in (9,17,33):
    x=np.linspace(0,1,n); y=x;X,Y=np.meshgrid(x,y)
    exact=.25*np.sin(np.pi*X)*np.cos(np.pi*Y)+.2*X+.1*Y
    ux=.25*np.pi*np.cos(np.pi*X)*np.cos(np.pi*Y)+.2
    source=exact*ux+1.4*.25*np.pi**2*np.sin(np.pi*X)*np.cos(np.pi*Y)
    u,h,r,ind=solve(.4,x,y,exact[:,0],exact[:,-1],np.full(n,.1),np.full(n,.1),source)
    refine.append({'nodes_per_axis':n,'max_error':float(np.max(abs(u-exact))),'iterations':len(h),'stencil_residual':ind})
orders=np.log2(np.array([v['max_error'] for v in refine[:-1]])/np.array([v['max_error'] for v in refine[1:]]))
assert np.all((orders>1.7)&(orders<2.3))
x=np.linspace(0,6,121); y=np.linspace(0,5,31);u=None;demo=[]
for eps in (.2,.1,.05):
    u,h,r,ind=solve(eps,x,y,np.ones(len(y)),-np.ones(len(y)),np.zeros(len(x)),np.zeros(len(x)),initial=u,tol=1e-8)
    demo.append({'epsilon':eps,'iterations':len(h),'matrix_residual':r,'stencil_residual':ind,'y_spread':float(np.max(np.ptp(u,axis=0)))})
result={'implementation':'Independent Python/SciPy cross-check; MATLAB not executed','quadratic':quad,'refinement':refine,'orders':orders.tolist(),'demo':demo}
print(json.dumps(result,indent=2))
