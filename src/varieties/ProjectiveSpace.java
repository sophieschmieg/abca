package varieties;

import java.util.ArrayList;
import java.util.Collections;
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
	public int getNumberOfElements() throws InfinityException {
		return this.field.getNumberOfElements() * (this.dimension + 1) - 1;
	}

	@Override
	public Iterable<ProjectivePoint<T>> getElements() throws InfinityException {
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
		/*
		List<Polynomial<T>> list = new ArrayList<Polynomial<T>>();
		int nonzero = p.getNonZero() - 1;
		List<T> coeff = new ArrayList<T>();
		for (int i = 0; i < this.dimension + 1; i++)
			coeff.add(this.field.zero());
		for (int i = 0; i < this.dimension + 1; i++) {
			if (i == nonzero)
				continue;
			coeff.set(i, this.field.one());
			coeff.set(nonzero, this.field.negative(p.getDehomogenisedCoord(i + 1, nonzero + 1)));
			list.add(this.ring.getLinear(coeff));
			coeff.set(i, this.field.zero());
		}
		return this.ring.getIdeal(list);*/
	}
	public PolynomialRing<T>.Ideal asHyperplaneIdeal(List<ProjectivePoint<T>> points) {
		if (points.isEmpty())
			return this.ring.getIdeal(Collections.singletonList(this.ring.one()));
		return this.asHyperplaneIdeal(points.subList(0, points.size() - 1)).intersect(this.asIdeal(points.get(0)));
	}
	public PolynomialRing<T>.Ideal asHyperplaneIdeal(ProjectivePoint<T> p1, ProjectivePoint<T> p2) {
		List<ProjectivePoint<T>> points = new ArrayList<ProjectivePoint<T>>();
		points.add(p1);
		points.add(p2);
		return this.asHyperplaneIdeal(points);
	}
}
