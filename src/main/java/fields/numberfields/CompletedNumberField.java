package fields.numberfields;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import fields.finitefields.FiniteField;
import fields.finitefields.FiniteField.FFE;
import fields.finitefields.PrimeField.PFE;
import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.helper.AbstractElement;
import fields.helper.AbstractFieldExtension;
import fields.helper.FieldEmbedding;
import fields.integers.FiniteRationalVectorSpace;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Integers.SmallestIntegerSolutionPreparation;
import fields.integers.LocalizedFractions;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.AlgebraicExtensionElement;
import fields.interfaces.DiscreteValuationField;
import fields.interfaces.DiscreteValuationRing.OkutsuType;
import fields.interfaces.FieldExtension;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.local.LocalRingExtension;
import fields.local.PAdicField;
import fields.local.PAdicField.PAdicNumber;
import fields.local.Value;
import fields.numberfields.CompletedNumberField.Ext;
import fields.numberfields.NumberField.NFE;
import fields.numberfields.NumberFieldIntegers.NumberFieldIdeal;
import fields.vectors.Matrix;
import fields.vectors.MatrixAlgebra;
import fields.vectors.Vector;
import util.Identity;
import util.MiscAlgorithms;

public class CompletedNumberField extends AbstractFieldExtension<PAdicNumber, Ext, CompletedNumberField>
		implements FieldExtension<PAdicNumber, Ext, CompletedNumberField>, DiscreteValuationField<Ext, FFE> {
	public static class Ext extends AbstractElement<Ext> implements AlgebraicExtensionElement<PAdicNumber, Ext> {
		private UnivariatePolynomial<Fraction> asPolynomial;
		private UnivariatePolynomial<PAdicNumber> asRoundedPolynomial;

		private Ext(UnivariatePolynomial<Fraction> asPolynomial,
				UnivariatePolynomial<PAdicNumber> asRoundedPolynomial) {
			this.asPolynomial = asPolynomial;
			this.asRoundedPolynomial = asRoundedPolynomial;
		}

		@Override
		public String toString() {
			return asRoundedPolynomial.toString(/*"α",*/ true);
		}

		@Override
		public int compareTo(Ext o) {
			return asPolynomial.compareTo(o.asPolynomial);
		}

		@Override
		public UnivariatePolynomial<PAdicNumber> asPolynomial() {
			return asRoundedPolynomial;
		}

	}

	private LocalizedNumberField localized;
	private OkutsuType<Fraction, PFE, PFE, FFE, FiniteField> type;
	private LocalRingExtension<PAdicNumber, PFE, Ext, CompletedNumberField, PFE, FFE, FiniteField> ringOfIntegers;
	private int accuracy;
	private PAdicField roundedBase;
	private LocalizedFractions exactBase;
	private Real uniformizerValue;
	private List<UnivariatePolynomial<Fraction>> integralBasis;
	private List<Integer> integralBasisValues;
	private Matrix<Fraction> toIntegralBaseChange;
	private MatrixAlgebra<Fraction> exactMatrixAlgebra;
	private Ext uniformizer;
	private OtherVersion<Ext, NFE, FFE, LocalizedNumberField> exact;
	private Map<Integer, CompletedNumberField> withAccuracy;
	private Map<Integer, SmallestIntegerSolutionPreparation> roundToIntegerPreparation;

	static CompletedNumberField getCompletedNumberField(LocalizedNumberField localized, int accuracy) {
		OkutsuType<Fraction, PFE, PFE, FFE, FiniteField> type = localized.type().clone();
		accuracy = MiscAlgorithms.DivRoundUp(accuracy, type.ramificationIndex()) * type.ramificationIndex();
		LocalizedFractions exactBase = Rationals.q().withValuation(localized.residueField().characteristic());
		type = exactBase.ringOfIntegers().singleFactorLifting(type, accuracy);
		PAdicField roundedBase = new PAdicField(exactBase.uniformizer().asInteger(),
				accuracy / type.ramificationIndex());
		UnivariatePolynomial<PAdicNumber> roundedPolynomial = roundedBase.getUnivariatePolynomialRing()
				.getEmbedding(type.representative(), roundedBase.fromRationalMap());
		return new CompletedNumberField(localized, accuracy, roundedBase, roundedPolynomial, type);
	}

	private CompletedNumberField(LocalizedNumberField localized, int accuracy, PAdicField roundedBase,
			UnivariatePolynomial<PAdicNumber> roundedPolynomial,
			OkutsuType<Fraction, PFE, PFE, FFE, FiniteField> type) {
		super(roundedPolynomial, roundedBase, localized.getVariableName());
		this.localized = localized;
		this.type = type;
		this.exactBase = Rationals.q().withValuation(roundedBase.getPrime().getValue());
		this.accuracy = accuracy;
		this.roundedBase = roundedBase;
		int degree = type.representative().degree();
		UnivariatePolynomialRing<Fraction> baseRing = exactBase.getUnivariatePolynomialRing();
		this.integralBasis = new ArrayList<>();
		this.integralBasisValues = new ArrayList<>();
		List<Vector<Fraction>> asVectors = new ArrayList<>();
		for (int i = 0; i < degree; i++) {
			UnivariatePolynomial<Fraction> integralPolynomial = type.divisorPolynomial(i);
			int value = type.valuation(integralPolynomial).value();
			int adjustedValue = MiscAlgorithms.DivRoundDown(value, type.ramificationIndex());
			integralPolynomial = baseRing.multiply(exactBase.power(exactBase.uniformizer(), -adjustedValue),
					integralPolynomial);
			integralBasis.add(integralPolynomial);
			integralBasisValues.add(value - adjustedValue * type.ramificationIndex());
			asVectors.add(baseRing.asVector(integralPolynomial, degree - 1));
		}
		Matrix<Fraction> fromIntegralBaseChange = Matrix.fromColumns(asVectors);
		this.exactMatrixAlgebra = new FiniteRationalVectorSpace(degree).matrixAlgebra();
		this.toIntegralBaseChange = exactMatrixAlgebra.inverse(fromIntegralBaseChange);
		List<Ext> integralBasisElements = new ArrayList<>();
		for (UnivariatePolynomial<Fraction> element : integralBasis) {
			integralBasisElements.add(fromExactPolynomial(element));
		}
		this.ringOfIntegers = new LocalRingExtension<>(this, roundedBase, type.reduction().extension(),
				integralBasisElements, Matrix.mapMatrix(roundedBase.fromRationalMap(), toIntegralBaseChange));
		this.uniformizer = fromExactPolynomial(type.lift(type.reduction().extension().one(), 1));
		this.withAccuracy = new TreeMap<>();
		this.roundToIntegerPreparation = new TreeMap<>();
	}

	@Override
	public Real value(Ext t) {
		Value v = valuation(t);
		if (v.isInfinite()) {
			return getReals().zero();
		}
		if (uniformizerValue == null) {
			uniformizerValue = getReals().power(roundedBase.value(roundedBase.uniformizer()),
					-type.ramificationIndex());
		}
		return getReals().power(uniformizerValue, -valuation(t).value());
	}

	@Override
	public Reals getReals() {
		return roundedBase.getReals();
	}

	@Override
	public boolean isComplete() {
		return true;
	}

	@Override
	public FactorizationResult<Polynomial<Ext>, Ext> factorization(UnivariatePolynomial<Ext> t) {
		return ringOfIntegers().factorization(t, true, accuracy);
	}

	@Override
	public Ext multiply(Ext t1, Ext t2) {
		return fromExactPolynomial(exactBase.getUnivariatePolynomialRing().multiply(t1.asPolynomial, t2.asPolynomial));
	}

	@Override
	public Ext inverse(Ext t) {
		return inverse(t, accuracy);
	}

	@Override
	public Ext inverse(Ext t, int accuracy) {
		if (accuracy > this.accuracy) {
			return withAccuracy(accuracy).inverse(t, accuracy);
		}
		Value valuation = valuation(t);
		if (valuation.isInfinite()) {
			throw new ArithmeticException("Division by zero!");
		}
		if (valuation.equals(Value.ZERO)) {
			UnivariatePolynomial<Fraction> inverse = exactBase.ringOfIntegers().invertInType(t.asPolynomial, type,
					accuracy);
			return fromSmallDegreeExactPolynomial(inverse);
		} else if (valuation.compareTo(Value.ZERO) < 0) {
			Ext uniformizerPower = fromSmallDegreeExactPolynomial(type.lift(residueField().one(), -valuation.value()));
			return multiply(inverse(multiply(t, uniformizerPower), accuracy), uniformizerPower);
		}
		Ext uniformizerPower = fromSmallDegreeExactPolynomial(type.lift(residueField().one(), -valuation.value()));
		return multiply(inverse(multiply(t, uniformizerPower), accuracy + valuation.value()), uniformizerPower);
	}

	@Override
	public Ext divide(Ext dividend, Ext divisor) {
		Value dividendValue = valuation(dividend);
		if (dividendValue.isInfinite()) {
			return zero();
		}
		Value divisorValue = valuation(divisor);
		if (divisorValue.isInfinite()) {
			throw new ArithmeticException("Division by zero!");
		}
		int neededAccuracy = accuracy + Math.abs(dividendValue.value()) + Math.abs(divisorValue.value());
		return multiply(dividend, inverse(divisor, neededAccuracy));
	}

	@Override
	public Ext add(Ext t1, Ext t2) {
		return fromSmallDegreeExactPolynomial(
				exactBase.getUnivariatePolynomialRing().add(t1.asPolynomial, t2.asPolynomial));
	}

	@Override
	public Ext negative(Ext t) {
		return fromSmallDegreeExactPolynomial(exactBase.getUnivariatePolynomialRing().negative(t.asPolynomial));
	}

	@Override
	public Ext negative(Ext t, int accuracy) {
		if (accuracy > this.accuracy) {
			withAccuracy(accuracy).negative(t, accuracy);
		}
		return round(negative(t), accuracy);
	}

	@Override
	public Value valuation(Ext t) {
		return type.valuation(t.asPolynomial);
	}

	@Override
	public Ext uniformizer() {
		return uniformizer;
	}

	@Override
	public FiniteField residueField() {
		return type.reduction().extension();
	}

	@Override
	public FFE reduceInteger(Ext t) {
		Value v = type.valuation(t.asPolynomial);
		if (v.compareTo(Value.ZERO) > 0) {
			return residueField().zero();
		}
		if (v.compareTo(Value.ZERO) < 0) {
			throw new ArithmeticException(t + " is not an integer");
		}
		return type.reduce(t.asPolynomial);
	}

	@Override
	public Ext upToUniformizerPower(Ext t) {
		Value v = type.valuation(t.asPolynomial);
		if (v.isInfinite()) {
			return zero();
		}
		return divide(t, power(uniformizer, v.value()));
	}

	@Override
	public Ext liftToInteger(FFE s) {
		return fromSmallDegreeExactPolynomial(type.lift(s, 0));
	}

	@Override
	public Ext round(Ext t, int accuracy) {
		return fromRoundedSmallDegreeExactPolynomial(roundPolynomial(t.asPolynomial, accuracy));
	}

	private UnivariatePolynomial<Fraction> roundPolynomial(UnivariatePolynomial<Fraction> t, int accuracy) {
		accuracy = Math.min(accuracy, this.accuracy);
		UnivariatePolynomialRing<Fraction> polynomialRing = exactBase.getUnivariatePolynomialRing();
		Vector<Fraction> asVector = polynomialRing.asVector(t, degree() - 1);
		Vector<Fraction> asIntegralVector = exactMatrixAlgebra.multiply(toIntegralBaseChange, asVector);
		UnivariatePolynomial<Fraction> result = polynomialRing.zero();
		for (int i = 0; i < degree(); i++) {
			Fraction coefficient = asIntegralVector.get(i + 1);
			Fraction rounded = exactBase.round(coefficient,
					MiscAlgorithms.DivRoundUp(accuracy - integralBasisValues.get(i), type.ramificationIndex()));
			result = polynomialRing.add(polynomialRing.multiply(rounded, integralBasis.get(i)), result);
		}
		return result;
	}

	@Override
	public int getAccuracy() {
		return accuracy;
	}

	@Override
	public CompletedNumberField withAccuracy(int accuracy) {
		if (!withAccuracy.containsKey(accuracy)) {
			withAccuracy.put(accuracy, getCompletedNumberField(localized, accuracy));
		}
		return withAccuracy.get(accuracy);
	}

	@Override
	public DiscreteValuationFieldExtension<Ext, FFE, PAdicNumber, Ext, FFE, CompletedNumberField> getUniqueExtension(
			UnivariatePolynomial<Ext> minimalPolynomial) {
		UnivariatePolynomial<NFE> exactPolynomial = localized.getUnivariatePolynomialRing()
				.getEmbedding(minimalPolynomial, exact().getRetraction());
		FieldEmbedding<Fraction, NFE, NumberField> embedding = localized.getNumberField()
				.getEmbeddedExtension(exactPolynomial);
		List<NumberFieldIdeal> ideals = embedding.getField().maximalOrder().idealsOver(localized.ideal(), embedding);
		if (ideals.size() != 1) {
			throw new ArithmeticException("Not irreducible");
		}
		CompletedNumberField result = getCompletedNumberField(
				embedding.getField().maximalOrder().localizeAndQuotient(ideals.get(0)), accuracy);
		return new DiscreteValuationFieldExtension<>(this, result, new MathMap<>() {
			@Override
			public Ext evaluate(Ext t) {
				NFE exactElement = exact().getRetraction().evaluate(t);
				NFE upperElement = embedding.getEmbedding(exactElement);
				return result.exact().getEmbedding().evaluate(upperElement);
			}
		}, new MathMap<>() {
			@Override
			public Vector<Ext> evaluate(Ext t) {
				NFE exactElement = result.exact().getRetraction().evaluate(t);
				Vector<NFE> asVector = embedding.asVector(exactElement);
				return Vector.mapVector(exact().getEmbedding(), asVector);
			}
		});
	}

	@Override
	public OtherVersion<Ext, NFE, FFE, LocalizedNumberField> exact() {
		if (exact == null) {
			NumberField field = localized.getNumberField();
			List<Vector<IntE>> asIntegerVectors = new ArrayList<>();
			Integers z = Integers.z();
			IntE prime = z.getInteger(localized.ideal().prime());
			IntE primePower = z.power(prime, MiscAlgorithms.DivRoundDown(accuracy, type.ramificationIndex()));
			MathMap<Fraction, IntE> asIntegerMap = new MathMap<>() {
				@Override
				public IntE evaluate(Fraction t) {
					return z.getInteger(t.getNumerator().getValue()
							.multiply(t.getDenominator().getValue().modInverse(primePower.getValue()))
							.mod(primePower.getValue()));
				}
			};
			for (NFE integral : field.maximalOrder().getModuleGenerators()) {
				Ext embedded = fromElement(integral);
				Vector<Fraction> asVector = exactMatrixAlgebra.multiply(toIntegralBaseChange,
						exactBase.getUnivariatePolynomialRing().asVector(embedded.asPolynomial, degree() - 1));
				asIntegerVectors.add(Vector.mapVector(asIntegerMap, asVector));
			}
			SmallestIntegerSolutionPreparation preparation = z.prepareSmallestIntegerSolution(asIntegerVectors,
					primePower);
			MathMap<Ext, NFE> retraction = new MathMap<>() {
				@Override
				public NFE evaluate(Ext t) {
					Vector<Fraction> asVector = exactMatrixAlgebra.multiply(toIntegralBaseChange,
							exactBase.getUnivariatePolynomialRing().asVector(t.asPolynomial, degree() - 1));
					Value valuation = Value.INFINITY;
					for (Fraction coeff : asVector.asList()) {
						valuation = valuation.min(exactBase.valuation(coeff));
					}
					if (valuation.compareTo(Value.ZERO) < 0) {
						List<Fraction> integers = new ArrayList<>();
						Fraction uniformizerPower = exactBase.power(exactBase.uniformizer(), -valuation.value());
						for (Fraction coeff : asVector.asList()) {
							integers.add(exactBase.multiply(coeff, uniformizerPower));
						}
						asVector = new Vector<>(integers);
					}
					Vector<IntE> asIntegerVector = Vector.mapVector(asIntegerMap, asVector);
					Vector<IntE> inIntegralBasis = z.smallestIntegerSolution(asIntegerVector, preparation);
					NFE result = field.maximalOrder().fromVector(inIntegralBasis);
					if (valuation.compareTo(Value.ZERO) < 0) {
						NFE uniformizerPowerInverse = field.power(field.getEmbedding(prime), valuation.value());
						result = field.multiply(uniformizerPowerInverse, result);
					}
					return result;
				}
			};
			this.exact = new OtherVersion<>(localized, retraction, fromElementMap());
		}
		return exact;
	}

	public Ext fromElement(NFE t) {
		return fromExactPolynomial(t.asPolynomial());
	}

	public MathMap<NFE, Ext> fromElementMap() {
		return new MathMap<>() {
			@Override
			public Ext evaluate(NFE t) {
				return fromElement(t);
			}
		};
	}

	public NFE roundToInteger(Ext t, int accuracy) {
		Integers z = Integers.z();
		accuracy = MiscAlgorithms.DivRoundUp(accuracy, type.ramificationIndex()) * type.ramificationIndex();
		t = round(t, accuracy);
		IntE primePower = z.power(localized.ideal().prime(), accuracy / type.ramificationIndex());
		MathMap<Fraction, IntE> asIntegerMap = new MathMap<>() {
			@Override
			public IntE evaluate(Fraction t) {
				return z.getInteger(t.getNumerator().getValue()
						.multiply(t.getDenominator().getValue().modInverse(primePower.getValue()))
						.mod(primePower.getValue()));
			}
		};
		NumberField field = localized.getNumberField();
		if (!roundToIntegerPreparation.containsKey(accuracy)) {
			List<Vector<IntE>> asIntegerVectors = new ArrayList<>();
			for (NFE integral : field.maximalOrder().getModuleGenerators()) {
				Ext embedded = fromElement(integral);
				Vector<Fraction> asVector = exactMatrixAlgebra.multiply(toIntegralBaseChange,
						exactBase.getUnivariatePolynomialRing().asVector(embedded.asPolynomial, degree() - 1));
				asIntegerVectors.add(Vector.mapVector(asIntegerMap, asVector));
			}
			roundToIntegerPreparation.put(accuracy, z.prepareSmallestIntegerSolution(asIntegerVectors, primePower));
		}
		Vector<Fraction> asVector = exactMatrixAlgebra.multiply(toIntegralBaseChange,
				exactBase.getUnivariatePolynomialRing().asVector(t.asPolynomial, degree() - 1));
		Value valuation = Value.INFINITY;
		for (Fraction coeff : asVector.asList()) {
			valuation = valuation.min(exactBase.valuation(coeff));
		}
		if (valuation.compareTo(Value.ZERO) < 0) {
			List<Fraction> integers = new ArrayList<>();
			Fraction uniformizerPower = exactBase.power(exactBase.uniformizer(), -valuation.value());
			for (Fraction coeff : asVector.asList()) {
				integers.add(exactBase.multiply(coeff, uniformizerPower));
			}
			asVector = new Vector<>(integers);
		}
		Vector<IntE> asIntegerVector = Vector.mapVector(asIntegerMap, asVector);
		Vector<IntE> inIntegralBasis = z.smallestIntegerSolution(asIntegerVector,
				roundToIntegerPreparation.get(accuracy));
		NFE result = field.maximalOrder().fromVector(inIntegralBasis);
		if (valuation.compareTo(Value.ZERO) < 0) {
			NFE uniformizerPowerInverse = field.power(field.getEmbedding(localized.ideal().prime()), valuation.value());
			result = field.multiply(uniformizerPowerInverse, result);
		}
		return result;
	}

	@Override
	public OtherVersion<Ext, Ext, FFE, CompletedNumberField> complete(int accuracy) {
		return new OtherVersion<>(withAccuracy(accuracy), new Identity<>(), new Identity<>());
	}

	@Override
	public Ext getRandomInteger() {
		List<PAdicNumber> random = new ArrayList<>();
		for (int i = 0; i < degree(); i++) {
			random.add(roundedBase.getRandomInteger());
		}
		return ringOfIntegers.fromVector(new Vector<>(random));
	}

	@Override
	public LocalRingExtension<PAdicNumber, PFE, Ext, CompletedNumberField, PFE, FFE, FiniteField> ringOfIntegers() {
		return ringOfIntegers;
	}

	@Override
	protected Ext fromSmallDegreePolynomial(UnivariatePolynomial<PAdicNumber> polynomial) {
		UnivariatePolynomial<Fraction> asExactPolynomial;
		if (exactBase == null) {
			OtherVersion<PAdicNumber, Fraction, PFE, LocalizedFractions> exact = ((PAdicField) getRing()).exact();
			asExactPolynomial = exact.getField().getUnivariatePolynomialRing().getEmbedding(polynomial,
					exact.getRetraction());
			return new Ext(asExactPolynomial, polynomial);
		}
		asExactPolynomial = exactBase.getUnivariatePolynomialRing().getEmbedding(polynomial,
				roundedBase.toRationalMap());
		return fromSmallDegreeExactPolynomial(asExactPolynomial);
	}

	private Ext fromSmallDegreeExactPolynomial(UnivariatePolynomial<Fraction> asPolynomial) {
		return fromRoundedSmallDegreeExactPolynomial(roundPolynomial(asPolynomial, accuracy));
	}

	private Ext fromRoundedSmallDegreeExactPolynomial(UnivariatePolynomial<Fraction> asPolynomial) {
		UnivariatePolynomial<PAdicNumber> asRoundedPolynomial = roundedBase.getUnivariatePolynomialRing()
				.getEmbedding(asPolynomial, roundedBase.fromRationalMap());
		return new Ext(asPolynomial, asRoundedPolynomial);
	}

	public Ext fromExactPolynomial(UnivariatePolynomial<Fraction> asPolynomial) {
		UnivariatePolynomialRing<Fraction> polynomialRing = exactBase.getUnivariatePolynomialRing();
		UnivariatePolynomial<Fraction> remainder = polynomialRing
				.toUnivariate(polynomialRing.remainder(asPolynomial, type.representative()));
		return fromSmallDegreeExactPolynomial(remainder);
	}

	@Override
	public Ext fromPolynomial(UnivariatePolynomial<PAdicNumber> polynomial) {
		if (exactBase == null) {
			UnivariatePolynomialRing<PAdicNumber> polynomialRing = getRing().getUnivariatePolynomialRing();
			return fromSmallDegreePolynomial(
					polynomialRing.toUnivariate(polynomialRing.remainder(polynomial, minimalPolynomial())));
		}
		return fromExactPolynomial(
				exactBase.getUnivariatePolynomialRing().getEmbedding(polynomial, roundedBase.toRationalMap()));
	}

	@Override
	public CompletedNumberField makeExtension(UnivariatePolynomial<PAdicNumber> minimalPolynomial) {
		UnivariatePolynomial<Fraction> exactMinimalPolynomial = exactBase.getUnivariatePolynomialRing()
				.getEmbedding(minimalPolynomial, roundedBase.toRationalMap());
		NumberField field = NumberField.getNumberField(exactMinimalPolynomial);
		List<NumberFieldIdeal> ideals = field.maximalOrder().idealsOver(residueField().characteristic());
		if (ideals.size() != 1) {
			throw new ArithmeticException("Not irreducible");
		}
		return getCompletedNumberField(field.maximalOrder().localizeAndQuotient(ideals.get(0)), accuracy);
	}

	@Override
	protected CompletedNumberField asExtensionType() {
		return this;
	}
}
