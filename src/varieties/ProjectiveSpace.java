package varieties;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import fields.Element;
import fields.Field;
import fields.InfinityException;
import fields.Polynomial;
import fields.PolynomialRing;

public class ProjectiveSpace<T extends Element> implements Variety<T> {
	private Field<T> field;
	private PolynomialRing<T> ring;
	private int dimension;

	public ProjectiveSpace(Field<T> field, int dimension) {
		this.field = field;
		this.dimension = dimension;
		this.ring = new PolynomialRing<T>(field, dimension + 1, Polynomial.GREVLEX);
	}

	public PolynomialRing<T> getRing() {
		return ring;
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
	public boolean hasRationalPoint(ProjectivePoint<T> p) {
		return p.getDim() == this.dimension;
	}

	public PolynomialRing<T>.Ideal asIdeal(ProjectivePoint<T> p) {
		return p.asIdeal(this.ring);
	}

	public PolynomialRing<T>.Ideal asHyperplaneIdeal(List<ProjectivePoint<T>> points) {
		if (points.isEmpty())
			return this.ring.getIdeal(Collections.singletonList(this.ring.one()));
		PolynomialRing<T>.Ideal intersection = this.asHyperplaneIdeal(points.subList(1, points.size()))
				.intersect(this.asIdeal(points.get(0)));
		List<Polynomial<T>> generators = new ArrayList<>();
		for (Polynomial<T> p : intersection.getBasis()) {
			if (p.getDegree() <= 1) {
				generators.add(p);
			}
		}
		return this.ring.getIdeal(generators);
	}

	public PolynomialRing<T>.Ideal asHyperplaneIdeal(ProjectivePoint<T> p1, ProjectivePoint<T> p2) {
		List<ProjectivePoint<T>> points = new ArrayList<ProjectivePoint<T>>();
		points.add(p1);
		points.add(p2);
		return this.asHyperplaneIdeal(points);
	}
}
