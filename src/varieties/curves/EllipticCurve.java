package varieties.curves;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.TreeMap;

import varieties.ProjectivePoint;
import varieties.ProjectiveSpace;
import varieties.ProjectiveVariety;
import varieties.curves.DivisorGroup.Divisor;
import fields.Element;
import fields.Field;
import fields.FunctionField;
import fields.Group;
import fields.InfinityException;
import fields.Polynomial;
import fields.PolynomialRing;
import fields.RationalFunction;
import fields.Ring;

public class EllipticCurve<T extends Element> extends ProjectiveVariety<T> implements SmoothCurve<T>, Group<ProjectivePoint<T>> {
	private Field<T> field;
	private T a;
	private T b;
	private ProjectivePoint<T> pointAtInfinity;
	private PolynomialRing<T> ring;
	private List<ProjectivePoint<T>> points;
	@SuppressWarnings("rawtypes")
	private static Map<Field, ProjectiveSpace> spaces = new HashMap<Field, ProjectiveSpace>();
        private Map<Integer, Polynomial<T>> divisionPolynomials;
	
	@SuppressWarnings("unchecked")
	private static<T extends Element> ProjectiveSpace<T> getSpace(Field<T> field) {
		if (!EllipticCurve.spaces.containsKey(field))
			EllipticCurve.spaces.put(field, new ProjectiveSpace<T>(field, 2));
		return EllipticCurve.spaces.get(field);
	}
	
	private static<T extends Element> PolynomialRing<T>.Ideal getIdeal(Field<T> field, T a, T b) {
		PolynomialRing<T> ring = EllipticCurve.getSpace(field).getRing();
		Polynomial<T> f = ring.getEmbedding(field.one(), new int[]{3, 0, 0});
		f = ring.add(f, ring.getEmbedding(field.negative(field.one()), new int[]{0, 2, 1}));
		f = ring.add(f, ring.getEmbedding(a, new int[]{1, 0, 2}));
		f = ring.add(f, ring.getEmbedding(b, new int[]{0, 0, 3}));
		return ring.getIdeal(Collections.singletonList(f));
	}
	@SuppressWarnings("unchecked")
	public EllipticCurve(Field<T> field, T a, T b) {
		super(EllipticCurve.getSpace(field), EllipticCurve.getIdeal(field, a, b));
		this.field = field;
		this.a = a;
		this.b = b;
		this.pointAtInfinity = new ProjectivePoint<T>(this.field, this.field.zero(), this.field.one(), this.field.zero());
		if (this.field.add(this.field.multiply(4, this.field.power(a, 3)), this.field.multiply(27, this.field.power(b, 2))).equals(this.field.zero()))
			throw new ArithmeticException("Singular curve");
		this.ring = this.getSpace().getRing();
	}
	/*@Override
	public boolean hasRationalPoint(ProjectivePoint<T> p) {
		if (p.getDim() != 2)
			return false;
		if (p.equals(this.pointAtInfinity))
			return true;
		if (p.getCoord(3).equals(this.field.zero()))
			return false;
		T leftside = this.field.power(p.getDehomogenisedCoord(2, 3), 2);
		T rightside =  this.field.power(p.getDehomogenisedCoord(1, 3), 3);
		rightside = this.field.add(rightside, this.field.multiply(this.a, p.getDehomogenisedCoord(1, 3)));
		rightside = this.field.add(rightside, this.b);
		if (leftside.equals(rightside))
			return true;
		return false;
	}*/

        public Polynomial<T> getDivisionPolynomial(int m) {
          if (this.divisionPolynomials == null) {
            this.divisionPolynomials = new TreeMap<>();
          }
          if (this.divisionPolynomials.containsKey(m)) {
            Polynomial<T> p = this.divisionPolynomials.get(m);
            if (m % 2 == 0) {
              return this.ring.multiply(2, this.ring.getVar(2), p);
            } else {
              return p;
            }
          }
          Polynomial<T> o = this.ring.one();
          Polynomial<T> x = this.ring.getVar(1);
          Polynomial<T> x2 = this.ring.getVar(1, 2);
          Polynomial<T> x3 = this.ring.getVar(1, 3);
          Polynomial<T> x4 = this.ring.getVar(1, 4);
          Polynomial<T> x6 = this.ring.getVar(1, 6);
          Polynomial<T> a = this.ring.getEmbedding(this.a);
          Polynomial<T> a2 = this.ring.multiply(a, a);
          Polynomial<T> a3 = this.ring.multiply(a2, a);
          Polynomial<T> b = this.ring.getEmbedding(this.b);
          Polynomial<T> b2 = this.ring.multiply(b, b);
          Polynomial<T> ab = this.ring.multiply(a, b);
          Polynomial<T> y2 = this.ring.add(x3, this.ring.multiply(a, x), b);
          Ring<Polynomial<T>> r = this.ring;
          if (m == 0) {
            this.divisionPolynomials.put(0, r.zero());
          } else if (m == 1) {
            this.divisionPolynomials.put(1, r.one());
          } else if (m == 2) {
            this.divisionPolynomials.put(2, r.one());
          } else if (m == 3) {
            Polynomial<T> p = r.multiply(3, x4);
            p = r.add(p, r.multiply(6, a, x2));
            p = r.add(p, r.multiply(12, b, x));
            p = r.add(p, r.multiply(-1, a2));
            this.divisionPolynomials.put(3, p);
          } else if (m == 4) {
            Polynomial<T> p = x6;
            p = r.add(p, r.multiply(5, a, x4));
            p = r.add(p, r.multiply(20, b, x3));
            p = r.add(p, r.multiply(-5, a2, x2));
            p = r.add(p, r.multiply(-4, ab, x));
            p = r.add(p, r.multiply(-8, b2));
            p = r.add(p, r.multiply(-1, a3));
            p = r.multiply(2, p);
            this.divisionPolynomials.put(4, p);
          } else if (m % 2 == 0) {
            int n = m / 2;
            for (int i = n - 2; i <= n + 2; i++) {
              this.getDivisionPolynomial(i);
            }
            Polynomial<T> psiNm2 = this.divisionPolynomials.get(n - 2);
            Polynomial<T> psiNm1 = this.divisionPolynomials.get(n - 1);
            Polynomial<T> psiN = this.divisionPolynomials.get(n);
            Polynomial<T> psiN1 = this.divisionPolynomials.get(n + 1);
            Polynomial<T> psiN2 = this.divisionPolynomials.get(n + 2);
            Polynomial<T> p = r.subtract(r.multiply(psiN2, r.power(psiNm1, 2)), r.multiply(psiNm2, r.power(psiN1, 2)));
            p = r.multiply(psiN, p);
            this.divisionPolynomials.put(m, p);
          } else if (m % 2 == 1) {
            int n = (m - 1) / 2;
            for (int i = n - 1; i <= n + 2; i++) {
              this.getDivisionPolynomial(i);
            }
            Polynomial<T> psiNm1 = this.divisionPolynomials.get(n - 1);
            Polynomial<T> psiN = this.divisionPolynomials.get(n);
            Polynomial<T> psiN1 = this.divisionPolynomials.get(n + 1);
            Polynomial<T> psiN2 = this.divisionPolynomials.get(n + 2);
            Polynomial<T> firstTerm = r.multiply(psiN2, r.power(psiN, 3));
            Polynomial<T> secondTerm = r.multiply(psiNm1, r.power(psiN1, 3));
            if (n % 2 == 0) {
              firstTerm = r.multiply(16, firstTerm, y2, y2);
            } else {
              secondTerm = r.multiply(16, secondTerm, y2, y2);
            }
            this.divisionPolynomials.put(m, r.subtract(firstTerm, secondTerm));
          }
          return this.getDivisionPolynomial(m);
        }

	public ProjectivePoint<T> neutral() {
		return this.pointAtInfinity;
	}
	@SuppressWarnings("unchecked")
	public ProjectivePoint<T> negative(ProjectivePoint<T> p) {
		return new ProjectivePoint<T>(this.field, p.getCoord(1), this.field.negative(p.getCoord(2)), p.getCoord(3));
	}
	public ProjectivePoint<T> add(ProjectivePoint<T> p, ProjectivePoint<T> q) {
		return this.negative(this.getThirdIntersection(p, q));
	}
	@Override
	public Field<T> getField() {
		return this.field;
	}
	@Override
	public int getEmbeddingDimension() {
		return 2;
	}
	@Override
	public boolean isProjective() {
		return true;
	}
	public RationalFunction<T> getRationalFunction(ProjectivePoint<T> zero, ProjectivePoint<T> pole1, ProjectivePoint<T> pole2) {
		ProjectivePoint<T> common = this.getThirdIntersection(pole1, pole2);
		Polynomial<T> numerator;
		Polynomial<T> denominator;
		if (common.equals(zero))
			numerator = this.getTangentSpace(zero).get(0);
		else
			numerator = this.getSpace().asHyperplaneIdeal(zero, common).getBasis().first();
		if (pole1.equals(pole2))
			denominator = this.getTangentSpace(pole1).get(0);
		else
			denominator = this.getSpace().asHyperplaneIdeal(pole1, pole2).getBasis().first();
		Set<ProjectivePoint<T>> zeroes = new TreeSet<ProjectivePoint<T>>();
		Set<ProjectivePoint<T>> poles = new TreeSet<ProjectivePoint<T>>();
		zeroes.add(zero);
		zeroes.add(this.getThirdIntersection(zero, common));
		poles.add(pole1);
		poles.add(pole2);
		return new RationalFunction<T>(this.field, numerator, denominator, this, this.ring);
	}
	@SuppressWarnings("unchecked")
	public ProjectivePoint<T> getThirdIntersection(ProjectivePoint<T> p, ProjectivePoint<T> q) {
		T s, r1, r2;
		if (p.equals(this.pointAtInfinity))
			return this.negative(q);
		if (q.equals(this.pointAtInfinity))
			return this.negative(p); 
		if (p.equals(this.negative(q)))
			return this.pointAtInfinity;
		T px = p.getDehomogenisedCoord(1, 3);
		T py = p.getDehomogenisedCoord(2, 3);
		T qx = q.getDehomogenisedCoord(1, 3);
		T qy = q.getDehomogenisedCoord(2, 3);
		if (px.equals(qx))
			s = this.field.divide(
					this.field.add(this.field.multiply(3, px, px), a),
					this.field.multiply(2, py));
		else
			s = this.field.divide(
					this.field.subtract(py, qy),
					this.field.subtract(px, qx));
		r1 = this.field.multiply(s, s);
		r1 = this.field.subtract(this.field.subtract(r1, px), qx);
		r2 = this.field.subtract(r1, px);
		r2 = this.field.multiply(s, r2);
		r2 = this.field.add(py, r2);
		return new ProjectivePoint<T>(this.field, r1, r2, this.field.one());
	}
	@Override
	public List<Polynomial<T>> getTangentSpace(ProjectivePoint<T> p) {
		List<T> line = new ArrayList<T>();
		List<Polynomial<T>> dl = getDifferentials();
		for (Polynomial<T> poly : dl)
			line.add(poly.evaluate(p.getCoords()));
		return Collections.singletonList(this.ring.getLinear(line));
	}
	@Override
	public List<Polynomial<T>> getCotangentSpace(ProjectivePoint<T> p) {
		List<Polynomial<T>> list = this.getDifferentials();
		List<Polynomial<T>> reslist = new ArrayList<Polynomial<T>>();;
		if (!this.pointAtInfinity.equals(p)) {
			reslist.add(list.get(1));
			reslist.add(this.ring.negative(list.get(0)));
			reslist.add(this.ring.zero());
		} else {
			reslist.add(list.get(2));
			reslist.add(this.ring.zero());
			reslist.add(this.ring.negative(list.get(0)));
		}
		return reslist;
	}
	private List<Polynomial<T>> getDifferentials() {
		List<Polynomial<T>> list = new ArrayList<Polynomial<T>>();
		Polynomial<T> x = this.ring.getVar(1);
		Polynomial<T> y = this.ring.getVar(2);
		Polynomial<T> z = this.ring.getVar(3);
		Polynomial<T> a = this.ring.getEmbedding(this.a);
		Polynomial<T> b = this.ring.getEmbedding(this.b);
		Ring<Polynomial<T>> r = this.ring;
		list.add(r.add(r.multiply(3, x, x), r.multiply(a, z, z)));
		list.add(r.multiply(-2, y, z));
		list.add(r.add(r.multiply(2, a, x, z), r.multiply(3, b, z, z), r.multiply(-1, y, y)));
		return list;
	}
	@SuppressWarnings("unchecked")
	@Override
	public ProjectivePoint<T> getRandomElement() {
		ProjectivePoint<T> p = null;
		do {
			T x = this.field.getRandomElement();
			T y = this.field.getRandomElement();
			T z = this.field.getRandomElement();
			if (x.equals(this.field.zero()) && y.equals(this.field.zero()) && z.equals(this.field.zero()))
				continue;
			p = new ProjectivePoint<T>(field, x, y, z);
		} while(!this.hasRationalPoint(p));
		return p;
	}
	@Override
	public boolean isFinite() {
		return this.field.isFinite();
	}
	@Override
	public int getNumberOfElements() {
		if (!this.field.isFinite())
			return -1;
                List<BigInteger> primes = new ArrayList<BigInteger>();
                BigInteger i = BigInteger.valueOf(3);
                BigInteger product = BigInteger.ONE;
                BigInteger q = this.field.getNumberOfElements();
                while (product.compareTo(BigInteger.valueOf(16).multiply(q)) < 0) {
                  primes.add(i);
                  product = product.multiply(i).multiply(i);
                  i = i.nextProbablePrime();
                }
                Polynomial<T> x = this.ring.getVar(1);
                Polynomial<T> x2 = this.ring.getVar(1, 2);
                Polynomial<T> x3 = this.ring.getVar(1, 3);
                Polynomial<T> y = this.ring.getVar(2);
                Polynomial<T> y2 = this.ring.getVar(2, 2);
                Polynomial<T> def = this.ring.add(x3, this.ring.multiply(this.a, x), this.ring.getEmbedding(b));
                def = this.ring.subtract(def, y2);
                for (BigInteger l : primes) {
                  PolynomialRing<T>.Ideal ideal = new PolynomialRing<T>.Ideal(def, this.getDivisionPolynomial(l.intValue()));
                  CoordinateRing<T> cr = new CoordinateRing<T>(this.ring, ideal);
                  CoordinateRingElement<T> xq = cr.power(cr.getEmbedding(x), q);
                  CoordinateRingElement<T> yq = cr.power(cr.getEmbedding(y), q);
                  CoordinateRingElement<T> xqq = cr.power(xq, q);
                  CoordinateRingElement<T> yqq = cr.power(yq, q);
                }

		int i = 0;
		Iterator<ProjectivePoint<T>> it;
		try {
			it = this.getElements().iterator();
		} catch (InfinityException e) {
			throw new RuntimeException(e);
		}
		while (it.hasNext()) {
			it.next();
			i++;
		}
		return i;
	}
	@SuppressWarnings("unchecked")
	@Override
	public Iterable<ProjectivePoint<T>> getElements() throws InfinityException {
		if (this.points != null)
			return Collections.unmodifiableList(this.points);
		this.points = new ArrayList<ProjectivePoint<T>>();
		this.points.add(this.pointAtInfinity);
		for (T x : this.field.getElements()) {
			for (T y : this.field.getElements()) {
				ProjectivePoint<T> p = new ProjectivePoint<T>(this.field, x, y, this.field.one());
				if (this.hasRationalPoint(p))
					this.points.add(p);
			}
		}
		return this.getElements();
	}
	@Override
	public ProjectivePoint<T> inverse(ProjectivePoint<T> t) {
		return this.negative(t);
	}
	@Override
	public ProjectivePoint<T> operate(ProjectivePoint<T> t1,
			ProjectivePoint<T> t2) {
		return this.add(t1, t2);
	}
        public ProjectivePoint<T> multiply(int n, ProjectivePoint<T> t) {
          ProjectivePoint<T> result = this.neutral();
          if (n < 0) {
            n = -n;
            t = this.inverse(t);
          }
          while (n > 0) {
            if ((n & 1) == 1) {
              result = this.operate(result, t);
            }
            n >>= 1;
            t = this.operate(t, t);
          }
          return result;
        }
	@Override
	public String toString() {
		return "Y^2Z = X^3 + " + this.a.toString() + "X + " + this.b.toString();
	}
	@Override
	public List<RationalFunction<T>> getRiemannRochSpace(Divisor<T> div) {
		if (div.getDegree() < 0)
			return Collections.emptyList();
		List<RationalFunction<T>> functions = new ArrayList<RationalFunction<T>>();
		FunctionField<T> ff = new FunctionField<T>(this.field, this, this.ring);
		List<ProjectivePoint<T>> zeroes = div.getPoles();
		List<ProjectivePoint<T>> poles = div.getZeroes();
		if (poles.size() == 0)
			return Collections.singletonList(ff.one());
		RationalFunction<T> f = ff.one();
		ProjectivePoint<T> firstpole = poles.get(0);
		ProjectivePoint<T> lastzero = null;
		int size = zeroes.size();
		if (zeroes.size() == poles.size()) {
			size--;
			lastzero = zeroes.get(size);
		}
		for (int i = 0; i < size; i++) {
			ProjectivePoint<T> zero = zeroes.get(i);
			ProjectivePoint<T> pole = poles.get(i + 1);
			f = ff.multiply(f, this.getRationalFunction(zero, firstpole, pole));
			firstpole = this.getThirdIntersection(firstpole, pole);
			firstpole = this.getThirdIntersection(firstpole, zero);
		}
		if (lastzero != null) {
			if (firstpole.equals(lastzero))
				return Collections.singletonList(f);
			else
				return Collections.emptyList();
		}
		functions.add(f);
		for (int i = size + 1; i < poles.size(); i++) {
			ProjectivePoint<T> pole = poles.get(i);
			ProjectivePoint<T> zero = this.getThirdIntersection(firstpole, pole);
			while (pole.equals(zero) || firstpole.equals(zero))
				zero = this.getRandomElement();
			functions.add(ff.multiply(f, this.getRationalFunction(zero, firstpole, pole)));
		}
		return functions;
	}
	@Override
	public boolean isPrincipal(Divisor<T> div) {
		if (div.getDegree() != 0)
			return false;
		return this.getRiemannRochSpace(div).size() == 1;
	}
}
