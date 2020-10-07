package fields.integers;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import fields.exceptions.InfinityException;
import fields.helper.AbstractElement;
import fields.helper.AbstractField;
import fields.helper.FieldOfFractions;
import fields.integers.Integers.IntE;
import fields.integers.Rationals.Fraction;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.numberfields.NumberField;
import fields.numberfields.NumberField.NFE;
import util.MiscAlgorithms;

public class Rationals extends AbstractField<Fraction> {
	private static Integers z = Integers.z();
	private static FieldOfFractions<IntE> q = new FieldOfFractions<>(z);
	private static Rationals rationals = new Rationals();
	
	public static class Fraction extends AbstractElement<Fraction>{
		private FieldOfFractions.Fraction<IntE> fraction;

		public Fraction(FieldOfFractions.Fraction<IntE> fraction) {
			this.fraction = fraction;
		}
		
		public FieldOfFractions.Fraction<IntE> fraction() {
			return fraction;
		}
		
		public IntE getNumerator() {
			return fraction.getNumerator();
		}
		

		public IntE getDenominator() {
			return fraction.getDenominator();
		}
		
		@Override
		public int compareTo(Fraction o) {
			return fraction.compareTo(o.fraction);
		}
		
		@Override
		public String toString() {
			return fraction.toString();
		}
	}
	
	private Rationals() {}
	
	public static Rationals q() {
		return rationals;
	}

	public List<Polynomial<Fraction>> factorization(Polynomial<Fraction> t) {
		if (t.numberOfVariables() != 1) {
			throw new RuntimeException("Multivariate factorization not implemented.");
		}
		PolynomialRing<Fraction> ring = this.getUnivariatePolynomialRing();
		NumberField asNumberField = new NumberField(ring.getVar(1));
		PolynomialRing<NFE> numberFieldRing = asNumberField.getUnivariatePolynomialRing();
		Polynomial<NFE> asNFEPoly = MiscAlgorithms.mapPolynomial(t, new MathMap<>() {
			@Override
			public NFE evaluate(Fraction n) {
				return asNumberField.getEmbedding(n);
			}}, numberFieldRing);
		List<Polynomial<NFE>> nfeResult = asNumberField.factorization(asNFEPoly);
		List<Polynomial<Fraction>> result = new ArrayList<>();
		for (Polynomial<NFE> nfeFactor : nfeResult) {
			result.add(MiscAlgorithms.mapPolynomial(nfeFactor, new MathMap<>() {
				@Override
				public Fraction evaluate(NFE n) {
					return asNumberField.asExtensionField().asPolynomial(n.asExtensionFieldElement()).coefficient(ring.getMonomial(new int[] {0}));
				}}, ring));
		}
		return result;
	}

	@Override
	public Fraction zero() {
		return new Fraction(q.zero());
	}

	@Override
	public Fraction one() {
		return new Fraction(q.one());
	}

	@Override
	public BigInteger characteristic() {
		return BigInteger.ZERO;
	}

	@Override
	public Fraction add(Fraction t1, Fraction t2) {
		return new Fraction(q.add(t1.fraction(), t2.fraction()));
	}

	@Override
	public Fraction negative(Fraction t) {
		return new Fraction(q.negative(t.fraction()));
	}

	@Override
	public Fraction multiply(Fraction t1, Fraction t2) {
		return new Fraction(q.multiply(t1.fraction(), t2.fraction()));
	}

	@Override
	public Fraction inverse(Fraction t) {
		return new Fraction(q.inverse(t.fraction()));
	}

	@Override
	public Fraction getRandomElement() {
		return new Fraction(q.getRandomElement());
	}

	@Override
	public boolean isFinite() {
		return false;
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		throw new InfinityException();
	}

	@Override
	public Iterator<Fraction> iterator() {
		return new Iterator<>() {
			Iterator<FieldOfFractions.Fraction<IntE>> it = q.iterator();

			@Override
			public boolean hasNext() {
				return it.hasNext();
			}

			@Override
			public Fraction next() {
				return new Fraction(it.next());
			}
		};
	}
	
	public Fraction getEmbedding(IntE t) {
		return new Fraction(q.getEmbedding(t));
	}
	
	public Fraction getFraction(IntE numerator, IntE denominator) {
		return new Fraction(q.getFraction(numerator, denominator));
	}
	
	public Fraction getEmbedding(BigInteger t) {
		return getEmbedding(z.getInteger(t));
	}
	
	public Fraction getInteger(BigInteger t) {
		return getEmbedding(t);
	}
	
	public boolean isInteger(Fraction t) {
		return t.fraction().getDenominator().equals(z.one());
	}
	
	@Override
	public String toString() {
		return "Q";
	}
}
