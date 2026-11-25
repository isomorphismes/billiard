"""Independent exact cellular oracle for the (pi/8, pi/2, 3pi/8) triangle.

No C kernel, generated reflection matrices, floating point, or third-party
packages. Run without arguments for the reference certificate and with an
ICK test executable to verify the emitted ledger and reject corrupted ledgers.
"""
import hashlib
import math
import subprocess
import sys
from collections import Counter, deque
from fractions import Fraction as Q
from pathlib import Path


def q(a=0, b=0):
    return Q(a), Q(b)


def add(a, b):
    return a[0] + b[0], a[1] + b[1]


def scale(a, n):
    return a[0]*n, a[1]*n


def mul(a, b):
    return a[0]*b[0]+2*a[1]*b[1], a[0]*b[1]+a[1]*b[0]


def plus(a, b):
    return tuple(add(x, y) for x, y in zip(a, b))


def times(a, n):
    return tuple(scale(x, n) for x in a)


ZERO = q(), q()
COS = (q(1), q(0, Q(1,2)), q(), q(0, -Q(1,2)),
       q(-1), q(0, -Q(1,2)), q(), q(0, Q(1,2)))
VERTICES = (ZERO, (q(1), q()), (q(1), q(-1,1)))


def image(chart, point):
    k, s = chart
    c, t = COS[k], COS[(k+6)%8]
    x, y = point
    return add(mul(c,x), scale(mul(t,y),-s)), add(mul(t,x), scale(mul(c,y),s))


def require(condition, reason):
    if not condition:
        raise ValueError(reason)


def model():
    # Right multiplication by the three parallel side reflections.
    charts = [(0,1)]
    seams = []
    for k, s in charts:
        row = []
        for j in (0,4,1):
            target = ((k+s*j)%8,-s)
            if target not in charts:
                charts.append(target)
            row.append(charts.index(target))
        seams.append(row)
    require(len(charts)==16, 'chart count')
    # Flood corner orbits (not the production union/find algorithm).
    classes, sizes = {}, []
    for f in range(16):
        for v in range(3):
            if (f,v) in classes:
                continue
            number = len(sizes)
            todo = deque([(f,v)])
            classes[f,v] = number
            count = 0
            while todo:
                face, vertex = todo.popleft()
                count += 1
                for side in (vertex,(vertex-1)%3):
                    corner = seams[face][side], vertex
                    if corner not in classes:
                        classes[corner] = number
                        todo.append(corner)
            sizes.append((v,count))
    require(Counter(sizes)==Counter({(0,16):1,(1,4):4,(2,16):1}), 'corner orbits')
    zeros = [i for i,(v,n) in enumerate(sizes) if v==2]
    edges, lookup = [], {}
    for f in range(16):
        for j,g in enumerate(seams[f]):
            require(seams[g][j]==f and charts[g][1]==-charts[f][1], 'seam involution')
            if f>g:
                continue
            a,b = VERTICES[j],VERTICES[(j+1)%3]
            h = image(charts[f],plus(b,times(a,-1)))
            require(h==image(charts[g],plus(b,times(a,-1))), 'seam displacement')
            lookup[f,j]=lookup[g,j]=len(edges)
            edges.append((f,j,g,classes[f,j],classes[f,(j+1)%3],h))
    require(len(edges)==24 and len(sizes)==6 and len(zeros)==1, 'cell counts')
    faces = []
    for f,(_,s) in enumerate(charts):
        chain = [0]*24
        for j in range(3):
            chain[lookup[f,j]] += s
        faces.append(chain)
    return edges,faces,zeros[0]


def boundary(edges, chain):
    result = [0]*6
    for (_,_,_,u,v,_),n in zip(edges,chain):
        result[u]-=n
        result[v]+=n
    return result


def holonomy(edges, chain):
    result = ZERO
    for edge,n in zip(edges,chain):
        result = plus(result,times(edge[-1],n))
    return result


def determinant(rows):
    """Exact elimination; only used on the small certificate matrix."""
    rows = [list(map(Q,row)) for row in rows]
    n = len(rows)
    require(all(len(row)==n for row in rows), 'non-square certificate')
    result = Q(1)
    for col in range(n):
        pivot = next((i for i in range(col,n) if rows[i][col]),None)
        if pivot is None:
            return Q(0)
        if pivot!=col:
            rows[pivot],rows[col]=rows[col],rows[pivot]
            result=-result
        p=rows[col][col]
        result*=p
        for i in range(col+1,n):
            factor=rows[i][col]/p
            for j in range(col+1,n):
                rows[i][j]-=factor*rows[col][j]
            rows[i][col]=0
    return result


def tree(edges, root, dual=False, excluded=()):
    paths={root:[0]*24}
    selected=set()
    todo=deque([root])
    while todo:
        v=todo.popleft()
        for i,(f,_,g,a,b,_) in enumerate(edges):
            if i in excluded:
                continue
            u,w=(f,g) if dual else (a,b)
            if v==u:
                neighbor,sgn=w,1
            elif v==w:
                neighbor,sgn=u,-1
            else:
                continue
            if neighbor in paths:
                continue
            selected.add(i)
            paths[neighbor]=paths[v].copy()
            paths[neighbor][i]+=sgn
            todo.append(neighbor)
    return selected,paths


def reference():
    edges,faces,root=model()
    primal,paths=tree(edges,root)
    dual,_=tree(edges,0,True,primal)
    periods={}
    for i,(_,_,_,u,v,_) in enumerate(edges):
        if i in primal or i in dual:
            continue
        chain=[a-b for a,b in zip(paths[u],paths[v])]
        chain[i]+=1
        periods[len(periods)]=(i,holonomy(edges,chain),chain)
    return edges,faces,root,primal,dual,periods


def check_certificate(edges,faces,root,primal,dual,periods):
    require(len(primal)==5 and len(dual)==15 and not primal&dual, 'tree sizes/disjointness')
    # Check actual selected edges connect the primal and dual graphs.
    require(len(tree(edges,root,excluded=set(range(24))-primal)[1])==6,'primal connectivity')
    require(len(tree(edges,0,True,set(range(24))-dual)[1])==16,'dual connectivity')
    require(len(periods)==4,'period count')
    require({record[0] for record in periods.values()}==set(range(24))-primal-dual,'witness set')
    require(all(sum(row[i] for row in faces)==0 for i in range(24)),'oriented face cancellation')
    for row in faces:
        require(boundary(edges,row)==[0]*6,'face boundary not closed')
        require(holonomy(edges,row)==ZERO,'face holonomy nonzero')
    for witness,value,row in periods.values():
        require(row[witness]==1 and all(n==0 for i,n in enumerate(row) if i!=witness and i not in primal), 'nonfundamental chain')
        require(boundary(edges,row)==[0]*6,'period not closed')
        require(holonomy(edges,row)==value,'incorrect exact period holonomy')
    # Integral coordinates on ker d1: retain non-primal coefficients.
    # Fifteen face boundaries (one omitted) plus four proposed cycles must
    # form a Z-basis of Z^19, not merely span over Q or R.
    non_tree=sorted(set(range(24))-primal)
    matrix=[[row[i] for i in non_tree] for row in faces[1:]+[r[2] for r in periods.values()]]
    det=determinant(matrix)
    require(abs(det)==1, f'not a saturated homology basis: det={det}')
    # Independently report rational coefficient rank and real rank.
    coeffs=[[c for scalar in value for c in scalar] for _,value,_ in periods.values()]
    qdet=determinant(coeffs)
    require(qdet!=0,'unexpected rational holonomy dependence for this fixture')
    vectors=[value for _,value,_ in periods.values()]
    require(any(add(mul(a[0],b[1]),scale(mul(a[1],b[0]),-1))!=q()
                for a in vectors for b in vectors),'unexpected real rank below two')
    return det,qdet


def decode(values):
    require(len(values)==6,'holonomy width')
    result=[]
    for a,b,d in (values[:3],values[3:]):
        require(d>0 and math.gcd(a,b,d)==1,'noncanonical exact number')
        result.append(q(Q(a,d),Q(b,d)))
    return tuple(result)


def verify_ledger(data):
    edges,faces,root=model()
    meta,actual,periods,chains={}, {}, {}, {}
    for line in data.decode('ascii').splitlines():
        name,*fields=line.split('\t')
        if name=='period_edge':
            i,f,j,g,u,v,p,d,*h=map(int,fields)
            require(i not in actual and 0<=i<24 and p in (0,1) and d in (0,1),'edge row')
            actual[i]=(f,j,g,u,v,decode(h),p,d)
        elif name=='period':
            i=int(fields[0]); require(fields[1]=='absolute' and i not in periods,'period row')
            periods[i]=(int(fields[2]),decode(list(map(int,fields[3:]))))
        elif name=='period_chain':
            i,e,n=map(int,fields)
            require(0<=i<4 and 0<=e<24 and n and (i,e) not in chains,'chain row')
            chains[i,e]=n
        elif name=='zero_order':
            require(name not in meta and fields==['0','2'],'zero order')
            meta[name]=fields
        else:
            require(name not in meta and len(fields)==1,'metadata row')
            meta[name]=fields[0]
    require(meta=={'period_radicand':'2','period_basis':'primal_tree_and_disjoint_dual_tree',
                  'genus':'2','positive_zeros':'1','absolute_periods':'4','relative_periods':'0',
                  'base_vertex':str(root),'zero_order':['0','2']},'field or topology metadata')
    require(set(actual)==set(range(24)) and set(periods)==set(range(4)),'missing records')
    for i,e in enumerate(edges):
        require(actual[i][:6]==e,f'edge {i}: wrong incidence or holonomy')
    primal={i for i in actual if actual[i][6]}
    dual={i for i in actual if actual[i][7]}
    complete={i:(*periods[i],[chains.get((i,e),0) for e in range(24)]) for i in range(4)}
    return check_certificate(edges,faces,root,primal,dual,complete)


def encode(point):
    result=[]
    for a,b in point:
        d=math.lcm(a.denominator,b.denominator)
        result.extend((int(a*d),int(b*d),d))
    return result


def reference_ledger(edges,root,primal,dual,periods):
    rows=[['period_radicand',2],['period_basis','primal_tree_and_disjoint_dual_tree'],
          ['genus',2],['positive_zeros',1],['absolute_periods',4],['relative_periods',0],
          ['base_vertex',root],['zero_order',0,2]]
    for i,(f,j,g,u,v,h) in enumerate(edges):
        rows.append(['period_edge',i,f,j,g,u,v,int(i in primal),int(i in dual),*encode(h)])
    for i,(w,h,chain) in periods.items():
        rows.append(['period',i,'absolute',w,*encode(h)])
        rows.extend(['period_chain',i,e,n] for e,n in enumerate(chain) if n)
    return ''.join('\t'.join(map(str,row))+'\n' for row in rows).encode('ascii')


def reject_mutants(data):
    rows=data.splitlines(keepends=True)
    period=next(row for row in rows if row.startswith(b'period\t'))
    fields=period.rstrip(b'\n').split(b'\t')
    fields[4]=str(int(fields[4])+int(fields[6])).encode('ascii')
    chain=next(row for row in rows if row.startswith(b'period_chain\t'))
    pieces=chain.rstrip(b'\n').split(b'\t')
    pieces[3]=str(2*int(pieces[3])).encode('ascii')
    mutants=(data.replace(b'period_radicand\t2\n',b'period_radicand\t3\n',1),
             data.replace(b'period_radicand\t2\n',b'',1),
             data.replace(period,b'\t'.join(fields)+b'\n',1),
             data.replace(chain,b'\t'.join(pieces)+b'\n',1))
    for i,mutant in enumerate(mutants):
        require(mutant!=data,'mutation did not apply')
        try:
            verify_ledger(mutant)
        except ValueError:
            continue
        raise ValueError(f'mutant {i} accepted')
    print('four ledger mutants rejected (wrong field, missing field, holonomy, chain)')


def main():
    edges,faces,root,primal,dual,periods=reference()
    det,qdet=check_certificate(edges,faces,root,primal,dual,periods)
    print(f'reference: V=6 E=24 F=16 g=2; det={det}; Q-coefficient det={qdet}; real holonomy rank=2')
    print('primal edges:',sorted(primal),'dual edges:',sorted(dual))
    for i,(w,h,chain) in periods.items():
        print(f'period {i}, witness {w}: {encode(h)}; chain {[(e,n) for e,n in enumerate(chain) if n]}')
    if len(sys.argv)==2:
        data=subprocess.check_output([str(Path(sys.argv[1]).resolve()),'--emit-eighth-periods'])
        label='ICK ledger input (compiler provenance must be checked separately)'
    elif len(sys.argv)==1:
        data=reference_ledger(edges,root,primal,dual,periods)
        label='reference-only ledger (no C execution)'
    else:
        raise SystemExit('usage: check_eighth_periods.py [test_periods_executable]')
    det,qdet=verify_ledger(data)
    print(f'{label}: certificate determinant={det}; rational coefficient determinant={qdet}')
    reject_mutants(data)
    print('ledger SHA256',hashlib.sha256(data).hexdigest())


if __name__=='__main__':
    main()
