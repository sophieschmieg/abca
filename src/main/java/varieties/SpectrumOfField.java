package varieties;

import java.math.BigInteger;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import fields.exceptions.InfinityException;
import fields.helper.AbstractElement;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.PolynomialRing;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.CoordinateRing;
import fields.polynomials.Monomial;
import varieties.SpectrumOfField.SingletonPoint;
import varieties.affine.AffineCover;
import varieties.affine.AffineMorphism;
import varieties.affine.AffinePoint;
import varieties.affine.AffineScheme;

public class SpectrumOfField<T extends Element<T>> extends AbstractScheme<T, SingletonPoint> {
	public static class SingletonPoint extends AbstractElement<SingletonPoint> {

		private SingletonPoint() {
		}

		@Override
		public int compareTo(SingletonPoint o) {
			return 0;
		}

		@Override
		public String toString() {
			return "O";
		}

	}

	public final static SingletonPoint POINT = new SingletonPoint();
	private Field<T> field;

	public SpectrumOfField(Field<T> field) {
		this.field = field;
	}

	@Override
	public String toString() {
		return "Spec " + field;
	}

	@Override
	public Exactness exactness() {
		return Exactness.EXACT;
	}

	@Override
	public SingletonPoint getRandomElement() {
		return POINT;
	}

	@Override
	public boolean isFinite() {
		return true;
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		return BigInteger.ONE;
	}

	@Override
	public Iterator<SingletonPoint> iterator() {
		return Collections.singleton(POINT).iterator();
	}

	@Override
	public Field<T> getField() {
		return field;
	}

	@Override
	public boolean hasRationalPoint(SingletonPoint p) {
		return true;
	}

	@Override
	public AffineCover<T> getAffineCover() {
		PolynomialRing<T> polynomials = AbstractPolynomialRing.getPolynomialRing(field, 0, Monomial.GREVLEX);
		AffineScheme<T> asAffineVariety = new AffineScheme<>(field, polynomials.getZeroIdeal().divideOut());
		return new AffineCover<>(Collections.singletonList(asAffineVariety),
				Collections.singletonList(Collections.singletonList(asAffineVariety.identityMorphism())));
	}

	@Override
	public List<Integer> affineCoverIndex(SingletonPoint p) {
		return Collections.singletonList(0);
	}

	@Override
	public AffinePoint<T> asAffinePoint(SingletonPoint p, int affineCoverIndex) {
		return new AffinePoint<>(field, Collections.emptyList());
	}

	@Override
	public SingletonPoint fromAffinePoint(AffinePoint<T> p, int affineCoverIndex) {
		return POINT;
	}

	@Override
	public Morphism<T, AffinePoint<T>, SingletonPoint> embedding(int affineCoverIndex) {
		return new Morphism<>() {

			@Override
			public SingletonPoint evaluate(AffinePoint<T> t) {
				return POINT;
			}

			@Override
			public AffineScheme<T> getDomain() {
				return getAffineCover().getCover().get(0);
			}

			@Override
			public SpectrumOfField<T> getRange() {
				return SpectrumOfField.this;
			}

			@Override
			public RestrictionResult<T> restrict(AffinePoint<T> preimage) {
				AffineMorphism<T> id = getAffineCover().getCover().get(0).identityMorphism();
				return new RestrictionResult<>(0, 0, id, id, id);
			}
		};
	}

	@Override
	public Morphism<T, SingletonPoint, SingletonPoint> pointAsMorphism(SingletonPoint p) {
		return new Morphism<>() {

			@Override
			public SingletonPoint evaluate(SingletonPoint t) {
				return POINT;
			}

			@Override
			public SpectrumOfField<T> getDomain() {
				return SpectrumOfField.this;
			}

			@Override
			public SpectrumOfField<T> getRange() {
				return SpectrumOfField.this;
			}

			@Override
			public RestrictionResult<T> restrict(SingletonPoint preimage) {
				AffineMorphism<T> id = getAffineCover().getCover().get(0).identityMorphism();
				return new RestrictionResult<>(0, 0, id, id, id);
			}
		};
	}

	@Override
	public List<Scheme<T, SingletonPoint>> irreducibleComponents() {
		return Collections.singletonList(this);
	}
}
