package varieties.curves;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import varieties.ProjectivePoint;
import varieties.curves.DivisorGroup.Divisor;
import fields.Element;
import fields.Field;
import fields.FunctionField;
import fields.InfinityException;
import fields.Polynomial;
import fields.PolynomialRing;
import fields.RationalFunction;

public class ProjectiveLine<T extends Element> implements SmoothCurve<T> {
	private Field<T> field;
	private ProjectivePoint<T> pointAtInfinity;
	private ProjectivePoint<T> pointAtZero;
	private PolynomialRing<T> ring;
	
	@SuppressWarnings("unchecked")
	public ProjectiveLine(Field<T> field) {
		this.field = field;
		this.ring = new PolynomialRing<T>(field, 2, Polynomial.GREVLEX);
		this.pointAtInfinity = new ProjectivePoint<T>(this.field, this.field.one(), this.field.zero());
		this.pointAtZero = new ProjectivePoint<T>(this.field, this.field.zero(), this.field.one());
	}

	@Override
	public Field<T> getField() {
		return this.field;
	}

	@Override
	public boolean hasRationalPoint(ProjectivePoint<T> p) {
		return p.getDim() == 1;
	}

	@SuppressWarnings("unchecked")
	@Override
	public ProjectivePoint<T> getRandomElement() {
		try {
			if (this.isFinite() && (int)(Math.random() * this.getNumberOfElements()) == 0)
				return this.pointAtInfinity;
		} catch (InfinityException e) {
			throw new RuntimeException(e);
		}
		return new ProjectivePoint<T>(this.field, this.field.getRandomElement(), this.field.one());
	}

	@Override
	public boolean isFinite() {
		return this.field.isFinite();
	}

	@Override
	public int getNumberOfElements() throws InfinityException {
		return this.field.getNumberOfElements() + 1;
	}

	@Override
	public Iterable<ProjectivePoint<T>> getElements() throws InfinityException {
		final Iterable<T> elements = field.getElements();
		return new Iterable<ProjectivePoint<T>>() {
			@Override
			public Iterator<ProjectivePoint<T>> iterator() {
				return new Iterator<ProjectivePoint<T>>() {
					private boolean first = true;
					private ProjectivePoint<T> nextpoint;
					private Iterator<T> xit = elements.iterator();
					
					@Override
					public boolean hasNext() {
						return first || nextpoint != null;
					}

					@SuppressWarnings("unchecked")
					@Override
					public ProjectivePoint<T> next() {
						ProjectivePoint<T> result = this.nextpoint;
						this.nextpoint = null;
						if (this.first) {
							result = pointAtInfinity;
							this.first = false;
						}
						if (xit.hasNext())
							this.nextpoint = new ProjectivePoint<T>(field, xit.next(), field.one());
						else
							this.nextpoint = null;
						return result;
					}

					@Override
					public void remove() {
						throw new UnsupportedOperationException();
					}
				};
			}
		};
	}

	@Override
	public int getEmbeddingDimension() {
		return 1;
	}

	@Override
	public boolean isProjective() {
		return true;
	}

	@Override
	public List<Polynomial<T>> getCotangentSpace(ProjectivePoint<T> p) {
		PolynomialRing<T> ring = new PolynomialRing<T>(field, 2, Polynomial.LEX);
		List<Polynomial<T>> coords = new ArrayList<Polynomial<T>>();
		for (int i = 1; i <= 2; i++)
			if (p.getNonZero() == i)
				coords.add(ring.zero());
			else
				coords.add(ring.one());
		return coords;
	}

	@Override
	public List<Polynomial<T>> getTangentSpace(ProjectivePoint<T> p) {
		return Collections.emptyList();
	}

	@SuppressWarnings("unchecked")
	@Override
	public List<RationalFunction<T>> getRiemannRochSpace(Divisor<T> div) {
		if (div.getDegree() < 0)
			return Collections.emptyList();
		List<RationalFunction<T>> functions = new ArrayList<RationalFunction<T>>();
		FunctionField<T> ff = new FunctionField<T>(this.field, this);
		List<ProjectivePoint<T>> zeroes = div.getPoles();
		List<ProjectivePoint<T>> poles = div.getZeroes();
		RationalFunction<T> f = ff.one();
		for (int i = 0; i < zeroes.size(); i++) {
			ProjectivePoint<T> zero = zeroes.get(i);
			ProjectivePoint<T> pole = poles.get(i);
			f = ff.multiply(f, new RationalFunction<T>(field, this.ring.getLinear(this.field.negative(zero.getCoord(2)), zero.getCoord(1)),
					this.ring.getLinear(this.field.negative(pole.getCoord(2)), pole.getCoord(1)), this));
		}
		functions.add(f);
		for (int i = zeroes.size(); i < poles.size(); i++) {
			ProjectivePoint<T> pole = poles.get(i);
			ProjectivePoint<T> zero = this.pointAtZero;
			if (pole.equals(this.pointAtZero))
				zero = this.pointAtInfinity;
			functions.add(ff.multiply(f, new RationalFunction<T>(field, 
					this.ring.getLinear(this.field.negative(zero.getCoord(2)), zero.getCoord(1)),
					this.ring.getLinear(this.field.negative(pole.getCoord(2)), pole.getCoord(1)),
					this)));
		}
		return functions;
	}

	@Override
	public boolean isPrincipal(Divisor<T> div) {
		return div.getDegree() == 0;
	}

}
