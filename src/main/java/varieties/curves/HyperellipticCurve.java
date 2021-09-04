/*package varieties.curves;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import varieties.ProjectivePoint;
import varieties.curves.DivisorGroup.Divisor;
import fields.Element;
import fields.Field;
import fields.InfinityException;
import fields.Polynomial;
import fields.PolynomialRing;
import fields.RationalFunction;
import fields.Ring;

public class HyperellipticCurve<T extends Element> implements SmoothCurve<T> {
	private Field<T> field;
	private List<T> lambda;
	private ProjectivePoint<T> pointAtInfinity;
	private Ring<Polynomial<T>> ring;
	private List<ProjectivePoint<T>> points;
	private Polynomial<T> definingpolynomial;
	
	@SuppressWarnings("unchecked")
	public HyperellipticCurve(Field<T> field, List<T> lambda) {
		this.field = field;
		this.lambda = new ArrayList<T>();
		this.lambda.addAll(lambda);
		this.pointAtInfinity = new ProjectivePoint<T>(this.field, this.field.zero(), this.field.one(), this.field.zero());
		if (this.lambda.size() % 2 != 1)
			throw new RuntimeException("Not supported");
		for (T l : this.lambda)
			if (l.equals(this.field.zero()) || l.equals(this.field.one()))
				throw new ArithmeticException("Singular curve");
		this.ring = new PolynomialRing<T>(this.field, 3, Polynomial.GREVLEX);
		Polynomial<T> x = Polynomial.getVar(this.field, 1, 3);
		Polynomial<T> y = Polynomial.getVar(this.field, 2, 3);
		Polynomial<T> z = Polynomial.getVar(this.field, 3, 3);
		this.definingpolynomial = this.ring.multiply(x, this.ring.substract(x, this.ring.one()));
		for (T l : this.lambda)
			this.definingpolynomial = this.ring.multiply(this.definingpolynomial, this.ring.substract(x, Polynomial.getEmbedding(this.field, l, 3)));
		this.definingpolynomial = this.ring.substract(this.definingpolynomial, this.ring.multiply(this.ring.power(y, 2), this.ring.power(z, this.lambda.size())));
	}
	@Override
	public boolean hasRationalPoint(ProjectivePoint<T> p) {
		if (p.getDim() != 2)
			return false;
		return this.definingpolynomial.evaluate(p.getCoords()).equals(this.field.zero());
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
	@Override
	public List<Polynomial<T>> getTangentSpace(ProjectivePoint<T> p) {
		List<T> line = new ArrayList<T>();
		List<Polynomial<T>> dl = getDifferentials();
		for (Polynomial<T> poly : dl)
			line.add(poly.evaluate(p.getCoords()));
		return Collections.singletonList(Polynomial.getLine(this.field, line));
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
		list.add(this.definingpolynomial.derivative(1));
		list.add(this.definingpolynomial.derivative(2));
		list.add(this.definingpolynomial.derivative(3));
		return list;
	}
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
	public String toString() {
		return this.definingpolynomial.toString();
	}
	@Override
	public List<RationalFunction<T>> getRiemannRochSpace(Divisor<T> div) {
		throw new RuntimeException("Not implemented");
	}
	@Override
	public boolean isPrincipal(Divisor<T> div) {
		if (div.getDegree() != 0)
			return false;
		return this.getRiemannRochSpace(div).size() == 1;
	}
}*/
