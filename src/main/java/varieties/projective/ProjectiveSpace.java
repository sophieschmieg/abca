package varieties.projective;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import fields.exceptions.InfinityException;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.Ideal;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.Monomial;

public class ProjectiveSpace<T extends Element<T>> extends AbstractProjectiveScheme<T> {
	private Field<T> field;
	private int dimension;

	public ProjectiveSpace(Field<T> field, int dimension) {
		super(new ProjectiveScheme<>(field,
				AbstractPolynomialRing.getPolynomialRing(field, dimension + 1, Monomial.GREVLEX),
				Collections.emptyList()));
		this.field = field;
		this.dimension = dimension;
	}

	@Override
	public ProjectivePoint<T> getRandomElement() {
		boolean notallnull = false;
		List<T> coords = new ArrayList<T>();
		do {
			coords.clear();
			for (int i = 0; i < this.dimension + 1; i++) {
				T t = this.field.getRandomElement();
				coords.add(t);
				if (!t.equals(this.field.zero()))
					notallnull = true;
			}
		} while (!notallnull);
		return new ProjectivePoint<T>(this.field, coords);
	}

	@Override
	public boolean isFinite() {
		return this.field.isFinite();
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		return this.field.getNumberOfElements().multiply(BigInteger.valueOf(this.dimension + 1))
				.subtract(BigInteger.ONE);
	}

	@Override
	public Iterator<ProjectivePoint<T>> iterator() {
		throw new UnsupportedOperationException();
	}

	@Override
	public Field<T> getField() {
		return this.field;
	}

	@Override
	public int dimension() {
		return dimension;
	}

	@Override
	public boolean hasRationalPoint(ProjectivePoint<T> p) {
		return p.getDim() == this.dimension;
	}

	public Ideal<Polynomial<T>> asIdeal(ProjectivePoint<T> p) {
		return p.asIdeal(asProjectiveVariety().homogenousPolynomialRing());
	}

	public Ideal<Polynomial<T>> asHyperplaneIdeal(List<ProjectivePoint<T>> points) {
		PolynomialRing<T> ring = asProjectiveVariety().homogenousPolynomialRing();
		if (points.isEmpty()) {
			return ring.getIdeal(Collections.singletonList(ring.one()));
		}
		Ideal<Polynomial<T>> intersection = ring.intersect(this.asHyperplaneIdeal(points.subList(1, points.size())),
				this.asIdeal(points.get(0)));
		List<Polynomial<T>> generators = new ArrayList<>();
		for (Polynomial<T> p : intersection.generators()) {
			if (p.degree() <= 1) {
				generators.add(p);
			}
		}
		return ring.getIdeal(generators);
	}

	public Ideal<Polynomial<T>> asHyperplaneIdeal(ProjectivePoint<T> p1, ProjectivePoint<T> p2) {
		List<ProjectivePoint<T>> points = new ArrayList<ProjectivePoint<T>>();
		points.add(p1);
		points.add(p2);
		return this.asHyperplaneIdeal(points);
	}

	@Override
	public List<ProjectiveSpace<T>> irreducibleComponents() {
		return Collections.singletonList(this);
	}
}
