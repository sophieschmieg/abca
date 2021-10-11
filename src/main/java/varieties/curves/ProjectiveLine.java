package varieties.curves;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import fields.exceptions.InfinityException;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.CoordinateRing;
import fields.polynomials.Monomial;
import varieties.FunctionField;
import varieties.RationalFunction;
import varieties.curves.DivisorGroup.Divisor;
import varieties.projective.AbstractProjectiveScheme;
import varieties.projective.ProjectivePoint;
import varieties.projective.GenericProjectiveScheme;

public class ProjectiveLine<T extends Element<T>> extends AbstractProjectiveScheme<T> implements SmoothCurve<T> {
	private Field<T> field;
	private ProjectivePoint<T> pointAtInfinity;
	private ProjectivePoint<T> pointAtZero;
	private PolynomialRing<T> ring;
	private CoordinateRing<T> coordRing;

	public ProjectiveLine(Field<T> field) {
		super(new GenericProjectiveScheme<>(field, AbstractPolynomialRing.getPolynomialRing(field, 2, Monomial.GREVLEX), Collections.emptyList()));
		this.field = field;
		this.ring = asGenericProjectiveScheme().homogenousPolynomialRing();
		PolynomialRing<T> r = field.getUnivariatePolynomialRing();
		this.coordRing = r.getZeroIdeal().divideOut();
		this.pointAtInfinity = new ProjectivePoint<T>(this.field, this.field.one(), this.field.zero());
		this.pointAtZero = new ProjectivePoint<T>(this.field, this.field.zero(), this.field.one());
		}

	@Override
	public Exactness exactness() {
		return field.exactness();
	}

	@Override
	public Field<T> getField() {
		return this.field;
	}

	@Override
	public boolean hasRationalPoint(ProjectivePoint<T> p) {
		return p.getDim() == 1;
	}

	@Override
	public ProjectivePoint<T> getRandomElement() {
		try {
			if (this.isFinite() && (int) (Math.random() * this.getNumberOfElements().intValue()) == 0)
				return this.pointAtInfinity;
		} catch (InfinityException e) {
			throw new RuntimeException(e);
		}
		return new ProjectivePoint<T>(this.field, this.field.getRandomElement(), this.field.one());
	}

//	@Override
//	public boolean hasSimplify() {
//		return true;
//	}
//
//	@Override
//	public List<CoordinateRingElement<T>> simplify(RationalFunction<T> t) {
//		PolynomialRing<T> r = coordRing.getPolynomialRing();
//		Polynomial<T> num = t.getNumerator().getElement();
//		List<CoordinateRingElement<T>> result = new ArrayList<>();
//		if (num.equals(r.zero())) {
//			result.add(coordRing.zero());
//			result.add(coordRing.one());
//			return result;
//		}
//		Polynomial<T> denom = t.getDenominator().getElement();
//		Polynomial<T> gcd = r.gcd(num, denom);
//		gcd = r.multiply(denom.leadingCoefficient(), gcd);
//		result.add(coordRing.getEmbedding(r.divide(num, gcd)));
//		result.add(coordRing.getEmbedding(r.divide(denom, gcd)));
//		return result;
//	}

	@Override
	public boolean isFinite() {
		return this.field.isFinite();
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		return this.field.getNumberOfElements().add(BigInteger.ONE);
	}

	@Override
	public Iterator<ProjectivePoint<T>> iterator() {
		return new Iterator<ProjectivePoint<T>>() {
			private boolean first = true;
			private ProjectivePoint<T> nextpoint;
			private Iterator<T> xit = field.iterator();

			@Override
			public boolean hasNext() {
				return first || nextpoint != null;
			}

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

	@Override
	public int getEmbeddingDimension() {
		return 1;
	}
	
	@Override
	public List<Polynomial<T>> getCotangentSpace(ProjectivePoint<T> p) {
		PolynomialRing<T> ring = AbstractPolynomialRing.getPolynomialRing(field, 2, Monomial.LEX);
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
		List<RationalFunction<T>> functions = new ArrayList<>();
		FunctionField<T> ff = new FunctionField<>(this);
		List<ProjectivePoint<T>> zeroes = div.getPoles();
		List<ProjectivePoint<T>> poles = div.getZeroes();
		RationalFunction<T> f = ff.one();
		for (int i = 0; i < zeroes.size(); i++) {
			ProjectivePoint<T> zero = zeroes.get(i);
			ProjectivePoint<T> pole = poles.get(i);
			f = ff.multiply(f,
					ff.getFunction(
							this.ring.getLinear(this.field.negative(zero.getCoord(2)), zero.getCoord(1)),
							this.ring.getLinear(this.field.negative(pole.getCoord(2)), pole.getCoord(1))));
		}
		functions.add(f);
		for (int i = zeroes.size(); i < poles.size(); i++) {
			ProjectivePoint<T> pole = poles.get(i);
			ProjectivePoint<T> zero = this.pointAtZero;
			if (pole.equals(this.pointAtZero))
				zero = this.pointAtInfinity;
			functions.add(ff.multiply(f,
					ff.getFunction(
							this.ring.getLinear(this.field.negative(zero.getCoord(2)), zero.getCoord(1)),
							this.ring.getLinear(this.field.negative(pole.getCoord(2)), pole.getCoord(1)))));
		}
		return functions;
	}

	@Override
	public boolean isPrincipal(Divisor<T> div) {
		return div.getDegree() == 0;
	}

	@Override
	public CoordinateRing<T> getCoordinateRing() {
		return this.coordRing;
	}

	@Override
	public List<ProjectiveLine<T>> irreducibleComponents() {
		return Collections.singletonList(this);
	}
	
}
