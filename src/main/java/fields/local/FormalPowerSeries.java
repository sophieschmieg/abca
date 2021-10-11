package fields.local;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

import fields.exceptions.InfinityException;
import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.helper.AbstractElement;
import fields.helper.AbstractField;
import fields.helper.FieldOfFractions.Fraction;
import fields.interfaces.AlgebraicExtensionElement;
import fields.interfaces.DiscreteValuationField;
import fields.interfaces.DiscreteValuationRing;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.FieldExtension;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.local.CompleteDVRExtension.Ext;
import fields.local.FormalPowerSeries.PowerSeries;
import fields.polynomials.CoordinateRing;
import fields.polynomials.LocalizedCoordinateRing;
import fields.polynomials.LocalizedCoordinateRing.LocalizedElement;
import fields.vectors.Matrix;
import fields.vectors.MatrixModule;
import fields.vectors.Vector;
import fields.vectors.pivot.PivotStrategy;
import fields.vectors.pivot.ValuationPivotStrategy;
import util.Identity;
import util.MiscAlgorithms;

public class FormalPowerSeries<T extends Element<T>> extends AbstractField<PowerSeries<T>>
		implements DiscreteValuationField<PowerSeries<T>, T> {
	private Field<T> field;
	private UnivariatePolynomialRing<T> ring;
	private int accuracy;
	Reals r = Reals.r(1024);
	private DiscreteValuationRing<PowerSeries<T>, T> localRing;
	private CoordinateRing<T> coordinateRing;
	private LocalizedCoordinateRing<T> localizedCoordinateRing;

	public static class PowerSeries<T extends Element<T>> extends AbstractElement<PowerSeries<T>> {
		private FormalPowerSeries<T> powerSeries;
		private UnivariatePolynomial<T> value;
		private int lowestPower;

		private PowerSeries(FormalPowerSeries<T> powerSeries, int lowestPower, UnivariatePolynomial<T> value) {
			this.powerSeries = powerSeries;
			if (lowestPower > powerSeries.accuracy) {
				value = powerSeries.ring.zero();
				lowestPower = 0;
			}
			value = powerSeries.ring.round(value, powerSeries.accuracy - lowestPower);
			if (value.equals(powerSeries.ring.zero())) {
				lowestPower = 0;
			}
			if (!value.order().isInfinite() && value.order().value() > 0) {
				lowestPower += value.order().value();
				value = powerSeries.ring.toUnivariate(
						powerSeries.ring.divide(value, powerSeries.ring.getVarPower(value.order().value())));
			}
			this.lowestPower = lowestPower;
			this.value = value;
		}

		public String toString() {
			StringBuilder buf = new StringBuilder().append("...");
			if (value.equals(powerSeries.ring.zero())) {
				return "0";
			}
			boolean first = true;
			for (int i = powerSeries.getAccuracy() - 1; i >= lowestPower; i--) {
				T digit = digit(i);
				if (digit.equals(powerSeries.field.zero())) {
					continue;
				}
				if (first) {
					first = false;
				} else {
					buf.append(" + ");
				}
				String coefficient = digit.toString();
				if (coefficient.contains(" ")) {
					coefficient = "(" + coefficient + ")";
				}
				if (i == 0) {
					buf.append(coefficient);
				} else {
					if (digit(i).equals(powerSeries.field.one())) {
						buf.append("ξ");
					} else {
						buf.append(coefficient + "*ξ");
					}
					if (i != 1) {
						buf.append("^" + (i));
					}
				}
			}
			return buf.toString();
		}

		public T digit(int position) {
			position -= lowestPower;
			if (position < 0) {
				return powerSeries.field.zero();
			}
			return value.univariateCoefficient(position);
		}

		@Override
		public int compareTo(PowerSeries<T> o) {
			PowerSeries<T> roundedThis = powerSeries.round(this, powerSeries.getAccuracy());
			PowerSeries<T> roundedOther = powerSeries.round(o, powerSeries.getAccuracy());
			boolean thisZero = false;
			if (roundedThis.value.equals(powerSeries.ring.zero())) {
				thisZero = true;
			}
			boolean otherZero = false;
			if (roundedOther.value.equals(powerSeries.ring.zero())) {
				otherZero = true;
			}
			if (thisZero && otherZero) {
				return 0;
			}
			if (thisZero) {
				return 1;
			}
			if (otherZero) {
				return -1;
			}
			if (roundedThis.lowestPower != roundedOther.lowestPower) {
				return roundedThis.lowestPower - roundedOther.lowestPower;
			}
			return roundedThis.value.compareTo(roundedOther.value);
		}
	}

	public FormalPowerSeries(Field<T> field, int accuracy) {
		accuracy = Math.max(accuracy, 3);
		this.accuracy = accuracy;
		this.field = field;
		this.ring = field.getUnivariatePolynomialRing();
		this.localRing = new LocalRingImplementation<>(this, field.toString() + "[[ξ]]");
		this.coordinateRing = ring.getZeroIdeal().divideOut();
		this.localizedCoordinateRing = new LocalizedCoordinateRing<>(field, coordinateRing,
				coordinateRing.getIdeal(Collections.singletonList(coordinateRing.getVar(1))));
	}

	@Override
	public Exactness exactness() {
		return Exactness.FIXED_POINT;
	}

	@Override
	public boolean isComplete() {
		return true;
	}

	@Override
	public DiscreteValuationFieldExtension<PowerSeries<T>, T, ?, ?, ?, ?> getUniqueExtension(
			UnivariatePolynomial<PowerSeries<T>> minimalPolynomial) {
		return getUniqueExtension(minimalPolynomial,
				residueField().getExtension(residueField().getUnivariatePolynomialRing().getVar()));
	}

	private <R extends Element<R>, RE extends AlgebraicExtensionElement<R, RE>, RFE extends FieldExtension<R, RE, RFE>> DiscreteValuationFieldExtension<PowerSeries<T>, T, PowerSeries<T>, Ext<PowerSeries<T>>, RE, CompleteDVRExtension<PowerSeries<T>, T, R, RE, RFE>> getUniqueExtension(
			UnivariatePolynomial<PowerSeries<T>> minimalPolynomial,
			Extension<T, R, RE, RFE> trivialReductionExtension) {
		CompleteDVRExtension<PowerSeries<T>, T, R, RE, RFE> extension = CompleteDVRExtension
				.getCompleteDVRExtension(minimalPolynomial, this, trivialReductionExtension);
		return new DiscreteValuationFieldExtension<>(this, extension, extension.getEmbeddingMap(),
				extension.asVectorMap());
	}

	@Override
	public Extension<PowerSeries<T>, ?, ?, ?> getExtension(UnivariatePolynomial<PowerSeries<T>> minimalPolynomial) {
		return getExtension(minimalPolynomial,
				residueField().getExtension(residueField().getUnivariatePolynomialRing().getVar()));
	}

	private <R extends Element<R>, RE extends AlgebraicExtensionElement<R, RE>, RFE extends FieldExtension<R, RE, RFE>> Extension<PowerSeries<T>, PowerSeries<T>, Ext<PowerSeries<T>>, CompleteDVRExtension<PowerSeries<T>, T, R, RE, RFE>> getExtension(
			UnivariatePolynomial<PowerSeries<T>> minimalPolynomial,
			Extension<T, R, RE, RFE> trivialReductionExtension) {
		CompleteDVRExtension<PowerSeries<T>, T, R, RE, RFE> extension = CompleteDVRExtension
				.getCompleteDVRExtension(minimalPolynomial, this, trivialReductionExtension);
		return new Extension<>(extension, this, extension.getEmbeddingMap(), extension.asVectorMap());
	}

//	private <R extends Element<R>, RE extends AlgebraicExtensionElement<R, RE>, RFE extends FieldExtension<R, RE, RFE>> Extension<PowerSeries<T>, PowerSeries<R>, Ext<PowerSeries<R>>, CompleteLocalFieldExtension<PowerSeries<R>, R, RE, RFE>> asExtension(
//			Extension<T, R, RE, RFE> reductionAsExtension) {
//		FormalPowerSeries<R> base = new FormalPowerSeries<>(reductionAsExtension.extension().getBaseField(), accuracy);
//		UnivariatePolynomialRing<PowerSeries<R>> polynomials = base.getUnivariatePolynomialRing();
//		CompleteLocalFieldExtension<PowerSeries<R>, R, RE, RFE> asExtension = new CompleteLocalFieldExtension<>(
//				polynomials.getEmbedding(reductionAsExtension.extension().minimalPolynomial(), new MathMap<>() {
//					@Override
//					public PowerSeries<R> evaluate(R t) {
//						return base.getEmbedding(t);
//					}
//				}), base, reductionAsExtension.extension());
//		return new Extension<>(asExtension, this, new MathMap<PowerSeries<T>, Ext<PowerSeries<R>>>() {
//
//			@Override
//			public Ext<PowerSeries<R>> evaluate(PowerSeries<T> t) {
//				List<PowerSeries<R>> coefficients = new ArrayList<>();
//				for (int i = 0; i < reductionAsExtension.extension().degree(); i++) {
//					coefficients.add(base.zero());
//				}
//				for (int j = 0; j <= t.value.degree(); j++) {
//					RE e = reductionAsExtension.embeddingMap().evaluate(t.value.univariateCoefficient(j));
//					Vector<R> asVector = reductionAsExtension.extension().asVector(e);
//					for (int i = 0; i < asVector.dimension(); i++) {
//						coefficients.set(i,
//								base.add(coefficients.get(i),
//										new PowerSeries<R>(base, t.lowestPower,
//												reductionAsExtension.extension().getBaseField()
//														.getUnivariatePolynomialRing()
//														.getEmbedding(asVector.get(i + 1), j))));
//					}
//				}
//				return asExtension.fromPolynomial(polynomials.getPolynomial(coefficients));
//			}
//		}, new MathMap<Ext<PowerSeries<R>>, Vector<PowerSeries<T>>>() {
//
//			@Override
//			public Vector<PowerSeries<T>> evaluate(Ext<PowerSeries<R>> t) {
//				Vector<PowerSeries<R>> asVector = asExtension.asVector(t);
//				PowerSeries<T> result = zero();
//				for (int i = 0; i < asVector.dimension(); i++) {
//					PowerSeries<R> powerSeries = asVector.get(i + 1);
//					for (int j = 0; j <= powerSeries.value.degree(); j++) {
//						R element = powerSeries.value.univariateCoefficient(j);
//						RE embeddedElement = reductionAsExtension.extension().fromPolynomial(reductionAsExtension
//								.extension().getBaseField().getUnivariatePolynomialRing().getEmbedding(element, i));
//						T asT = reductionAsExtension.asVectorMap().evaluate(embeddedElement).get(1);
//						PowerSeries<T> asPowerSeries = new PowerSeries<>(FormalPowerSeries.this,
//								powerSeries.lowestPower, ring.getEmbedding(asT, j));
//						result = add(result, asPowerSeries);
//					}
//				}
//				return new Vector<>(Collections.singletonList(result));
//			}
//		});
//	}
//
//	private <R extends Element<R>, RE extends AlgebraicExtensionElement<R, RE>, RFE extends FieldExtension<R, RE, RFE>> Extension<PowerSeries<T>, PowerSeries<R>, Ext<PowerSeries<R>>, CompleteLocalFieldExtension<PowerSeries<R>, R, RE, RFE>> getExtension(
//			UnivariatePolynomial<PowerSeries<T>> minimalPolynomial,
//			Extension<PowerSeries<T>, PowerSeries<R>, Ext<PowerSeries<R>>, CompleteLocalFieldExtension<PowerSeries<R>, R, RE, RFE>> asTrivialExtension) {
//		UnivariatePolynomialRing<Ext<PowerSeries<R>>> basePolynomials = asTrivialExtension.extension()
//				.getUnivariatePolynomialRing();
//		Extension<Ext<PowerSeries<R>>, PowerSeries<R>, Ext<PowerSeries<R>>, CompleteLocalFieldExtension<PowerSeries<R>, R, RE, RFE>> extension = asTrivialExtension
//				.extension()
//				.getExtension(basePolynomials.getEmbedding(minimalPolynomial, asTrivialExtension.embeddingMap()));
//		return new Extension<>(extension.extension(), this,
//				new ConCatMap<>(asTrivialExtension.embeddingMap(), extension.embeddingMap()),
//				new ConCatMap<Ext<PowerSeries<R>>, Vector<Ext<PowerSeries<R>>>, Vector<PowerSeries<T>>>(
//						extension.asVectorMap(), new MathMap<>() {
//							@Override
//							public Vector<PowerSeries<T>> evaluate(Vector<Ext<PowerSeries<R>>> t) {
//								return asTrivialExtension.asVectorMap().evaluate(t.get(1));
//							}
//						}));
//
//	}

	public int getAccuracy() {
		return accuracy;
	}

	public FormalPowerSeries<T> withAccuracy(int accuracy) {
		return new FormalPowerSeries<>(field, accuracy);
	}

	@Override
	public OtherVersion<PowerSeries<T>, LocalizedElement<T>, T, LocalizedCoordinateRing<T>> exact() {
		return new OtherVersion<>(localizedCoordinateRing, toRationalFunctionMap(), fromRationalFunctionMap());
	}

	@Override
	public OtherVersion<PowerSeries<T>, PowerSeries<T>, T, FormalPowerSeries<T>> complete(int accuracy) {
		return new OtherVersion<>(withAccuracy(accuracy), new Identity<>(), new Identity<>());
	}

	@Override
	public String toString() {
		return field.toString() + "((ξ))";
	}

	@Override
	public Real value(PowerSeries<T> t) {
		if (t.equals(zero())) {
			return r.zero();
		}
		return r.exp(r.getInteger(-valuation(t).value()));
	}

	@Override
	public PowerSeries<T> zero() {
		return getEmbedding(ring.zero());
	}

	@Override
	public PowerSeries<T> one() {
		return getEmbedding(ring.one());
	}

	@Override
	public BigInteger characteristic() {
		return field.characteristic();
	}

	@Override
	public PowerSeries<T> add(PowerSeries<T> t1, PowerSeries<T> t2) {
		int lowestPower = Math.min(t1.lowestPower, t2.lowestPower);
		UnivariatePolynomial<T> p1 = ring.multiply(t1.value, ring.getVarPower(t1.lowestPower - lowestPower));
		UnivariatePolynomial<T> p2 = ring.multiply(t2.value, ring.getVarPower(t2.lowestPower - lowestPower));
		return getElement(ring.add(p1, p2), lowestPower);
	}

	@Override
	public PowerSeries<T> negative(PowerSeries<T> t) {
		return getElement(ring.negative(t.value), t.lowestPower);
	}

	@Override
	public PowerSeries<T> negative(PowerSeries<T> t, int accuracy) {
		return negative(t);
	}

	@Override
	public PowerSeries<T> multiply(PowerSeries<T> t1, PowerSeries<T> t2) {
		return getElement(ring.multiply(t1.value, t2.value), t1.lowestPower + t2.lowestPower);
	}

	@Override
	public PowerSeries<T> inverse(PowerSeries<T> t) {
		return divide(one(), t);
	}

	@Override
	public PowerSeries<T> inverse(PowerSeries<T> t, int accuracy) {
		if (accuracy < this.accuracy) {
			return round(inverse(t), accuracy);
		} else if (accuracy == this.accuracy) {
			return inverse(t);
		}
		return withAccuracy(accuracy).inverse(t);
	}

	@Override
	public PowerSeries<T> divide(PowerSeries<T> dividend, PowerSeries<T> divisor) {
		int lowestPower = dividend.lowestPower - divisor.lowestPower;
		List<T> coefficients = new ArrayList<>();
		for (int i = 0; i <= accuracy - lowestPower; i++) {
			T value = dividend.value.univariateCoefficient(i);
			for (int j = 0; j < i; j++) {
				value = field.subtract(value,
						field.multiply(coefficients.get(j), divisor.value.univariateCoefficient(i - j)));
			}
			coefficients.add(field.divide(value, divisor.value.univariateCoefficient(0)));
		}
		return getElement(ring.getPolynomial(coefficients), lowestPower);
	}

	@Override
	public PivotStrategy<PowerSeries<T>> preferredPivotStrategy() {
		return new ValuationPivotStrategy<>(this.valuation());
	}

	@Override
	public PowerSeries<T> getRandomElement() {
		int lowestPower = new Random().nextInt(getAccuracy());
		lowestPower -= getAccuracy() / 2;
		return getElement(ring.toUnivariate(ring.getRandomElement()), lowestPower);
	}

	@Override
	public PowerSeries<T> getRandomInteger() {
		return getElement(ring.toUnivariate(ring.getRandomElement()), 0);
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
	public Iterator<PowerSeries<T>> iterator() {
		throw new InfinityException();
	}

	@Override
	public boolean isInteger(PowerSeries<T> t) {
		return t.lowestPower >= 0;
	}

	@Override
	public Value valuation(PowerSeries<T> t) {
		if (t.value.equals(ring.zero())) {
			return Value.INFINITY;
		}
		return new Value(t.lowestPower);
	}

	public PowerSeries<T> getEmbedding(T value) {
		return getEmbedding(ring.getEmbedding(value));
	}

	public PowerSeries<T> getEmbedding(UnivariatePolynomial<T> value) {
		return getElement(value, 0);
	}

	public PowerSeries<T> getEmbedding(Fraction<Polynomial<T>> value) {
		return divide(getEmbedding(ring.toUnivariate(value.getNumerator())),
				getEmbedding(ring.toUnivariate(value.getDenominator())));
	}

	public PowerSeries<T> getElement(UnivariatePolynomial<T> value, int lowestPower) {
		return new PowerSeries<T>(this, lowestPower, value);
	}

	@Override
	public PowerSeries<T> uniformizer() {
		return getEmbedding(ring.getVar());
	}

	@Override
	public Field<T> residueField() {
		return field;
	}

	@Override
	public T reduceInteger(PowerSeries<T> t) {
		if (t.lowestPower < 0) {
			throw new ArithmeticException("not integer");
		}
		if (t.lowestPower > 0) {
			return field.zero();
		}
		return t.value.univariateCoefficient(0);
	}

	@Override
	public PowerSeries<T> upToUniformizerPower(PowerSeries<T> t) {
		return getEmbedding(t.value);
	}

	@Override
	public PowerSeries<T> liftToInteger(T s) {
		return getEmbedding(ring.getEmbedding(s));
	}

	@Override
	public PowerSeries<T> round(PowerSeries<T> t, int accuracy) {
		if (t.lowestPower > accuracy) {
			return zero();
		}
		return getElement(ring.round(t.value, accuracy - t.lowestPower), t.lowestPower);
	}

	public UnivariatePolynomial<T> roundToPolynomial(PowerSeries<T> t, int accuracy) {
		if (t.lowestPower < 0) {
			throw new ArithmeticException("Not near a polynomial");
		}
		t = round(t, accuracy);
		return ring.multiply(t.value, ring.getVarPower(t.lowestPower));
	}

	public UnivariatePolynomial<T> toPolynomial(PowerSeries<T> t) {
		return roundToPolynomial(t, getAccuracy());
	}

	public MathMap<PowerSeries<T>, Polynomial<T>> toPolynomialMap() {
		return new MathMap<>() {

			@Override
			public Polynomial<T> evaluate(PowerSeries<T> t) {
				return toPolynomial(t);
			}
		};
	}

	public PowerSeries<T> fromPolynomial(UnivariatePolynomial<T> t) {
		return getEmbedding(t);
	}

	public MathMap<Polynomial<T>, PowerSeries<T>> fromPolynomialMap() {
		return new MathMap<>() {

			@Override
			public PowerSeries<T> evaluate(Polynomial<T> t) {
				return fromPolynomial(ring.toUnivariate(t));
			}
		};
	}

	public LocalizedElement<T> roundToRationalFunction(PowerSeries<T> t) {
		int lowestPower = t.lowestPower;
		int balancedLength = MiscAlgorithms.DivRoundUp(accuracy - lowestPower, 2);
		int minDigit = accuracy - balancedLength * 2;
		List<List<T>> equations = new ArrayList<>();
		List<T> rhs = new ArrayList<>();
		for (int i = 0; i < balancedLength; i++) {
			int digit = i + minDigit;
			rhs.add(t.digit(digit));
			List<T> row = new ArrayList<>();
			for (int j = 0; j < i; j++) {
				row.add(field.zero());
			}
			row.add(field.one());
			for (int j = i + 1; j < balancedLength; j++) {
				row.add(field.zero());
			}
			for (int j = 1; j <= i; j++) {
				row.add(field.negative(t.digit(digit - j)));
			}
			for (int j = i; j < balancedLength; j++) {
				row.add(field.zero());
			}
			equations.add(row);
		}
		for (int i = 0; i < balancedLength; i++) {
			int digit = i + minDigit + balancedLength;
			rhs.add(t.digit(digit));
			List<T> row = new ArrayList<>();
			for (int j = 0; j < balancedLength; j++) {
				row.add(field.zero());
			}
			for (int j = 1; j <= balancedLength; j++) {
				row.add(field.negative(t.digit(digit - j)));
			}
			equations.add(row);
		}
		Matrix<T> matrix = new Matrix<>(equations);
		Vector<T> b = new Vector<>(rhs);
		Vector<T> x = new MatrixModule<>(field, 2 * balancedLength, 2 * balancedLength).solve(matrix, b);
		Polynomial<T> numerator = ring.getPolynomial(x.asList().subList(0, balancedLength));
		Polynomial<T> denominator = ring.add(ring.one(), ring
				.multiply(ring.getPolynomial(x.asList().subList(balancedLength, 2 * balancedLength)), ring.getVar()));
		numerator = ring.multiply(numerator, ring.getVarPower(Math.max(minDigit, 0)));
		denominator = ring.multiply(denominator, ring.getVarPower(Math.max(-minDigit, 0)));
		return localizedCoordinateRing.getEmbedding(numerator, denominator);
	}

	public MathMap<PowerSeries<T>, LocalizedElement<T>> toRationalFunctionMap() {
		return new MathMap<>() {

			@Override
			public LocalizedElement<T> evaluate(PowerSeries<T> t) {
				return roundToRationalFunction(t);
			}
		};
	}

	public PowerSeries<T> fromRationalFunction(LocalizedElement<T> t) {
		UnivariatePolynomial<T> numerator = ring.toUnivariate(t.getNumerator().getElement());
		UnivariatePolynomial<T> denominator = ring.toUnivariate(t.getDenominator().getElement());
		return divide(getEmbedding(numerator), getEmbedding(denominator));
	}

	public MathMap<LocalizedElement<T>, PowerSeries<T>> fromRationalFunctionMap() {
		return new MathMap<>() {
			@Override
			public PowerSeries<T> evaluate(LocalizedElement<T> t) {
				return fromRationalFunction(t);
			}
		};
	}

	@Override
	public DiscreteValuationRing<PowerSeries<T>, T> ringOfIntegers() {
		return localRing;
	}

	@Override
	public PowerSeries<T> characteristicRoot(PowerSeries<T> t, int power) {
		return getElement(ring.toUnivariate(ring.characteristicRoot(t.value, power)),
				BigInteger.valueOf(t.lowestPower).divide(characteristic().pow(power)).intValueExact());
	}

	@Override
	public boolean hasCharacteristicRoot(PowerSeries<T> t, int power) {
		return ring.hasCharacteristicRoot(t.value)
				&& BigInteger.valueOf(t.lowestPower).mod(characteristic().pow(power)).equals(BigInteger.ZERO);
	}

	@Override
	public boolean isIrreducible(UnivariatePolynomial<PowerSeries<T>> t) {
		return ringOfIntegers().isIrreducible(t);
	}

	@Override
	public FactorizationResult<Polynomial<PowerSeries<T>>, PowerSeries<T>> factorization(
			UnivariatePolynomial<PowerSeries<T>> t) {
		return ringOfIntegers().factorization(t, true, getAccuracy());
	}
}
