package fields.local;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.helper.AbstractElement;
import fields.helper.AbstractFieldExtension;
import fields.interfaces.AlgebraicExtensionElement;
import fields.interfaces.DiscreteValuationField;
import fields.interfaces.DiscreteValuationField.ExtensionOfDiscreteValuationField;
import fields.interfaces.DiscreteValuationField.OtherVersion;
import fields.interfaces.DiscreteValuationFieldExtension;
import fields.interfaces.DiscreteValuationRing.OkutsuType;
import fields.interfaces.DiscreteValuationRing.TheMontesResult;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.FieldExtension;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.local.CompleteDVRExtension.Ext;
import fields.vectors.Matrix;
import fields.vectors.Vector;
import util.Identity;
import util.MiscAlgorithms;

public class CompleteDVRExtension<B extends Element<B>, S extends Element<S>, R extends Element<R>, RE extends AlgebraicExtensionElement<R, RE>, RFE extends FieldExtension<R, RE, RFE>>
		extends AbstractFieldExtension<B, Ext<B>, CompleteDVRExtension<B, S, R, RE, RFE>>
		implements DiscreteValuationFieldExtension<B, S, Ext<B>, CompleteDVRExtension<B, S, R, RE, RFE>, R, RE, RFE> {
	public static class Ext<B extends Element<B>> extends AbstractElement<Ext<B>>
			implements AlgebraicExtensionElement<B, Ext<B>> {
		private UnivariatePolynomial<B> asPolynomial;
		private CompleteDVRExtension<B, ?, ?, ?, ?> extension;

		private Ext(CompleteDVRExtension<B, ?, ?, ?, ?> extension, UnivariatePolynomial<B> asPolynomial) {
			this.extension = extension;
			this.asPolynomial = asPolynomial;
		}

		@Override
		public int compareTo(Ext<B> o) {
			Ext<B> roundedThis = extension.round(this, extension.getAccuracy());
			Ext<B> roundedOther = extension.round(o, extension.getAccuracy());
			return roundedThis.asPolynomial.compareTo(roundedOther.asPolynomial);
		}

		@Override
		public UnivariatePolynomial<B> asPolynomial() {
			return asPolynomial;
		}

		@Override
		public String toString() {
			DiscreteValuationField<B, ?> lowAccuracyBase = extension.getBaseField()
					.withAccuracy(MiscAlgorithms.DivRoundDown(extension.getAccuracy(), extension.ramificationIndex()));
			UnivariatePolynomialRing<B> univariatePolynomialRing = lowAccuracyBase.getUnivariatePolynomialRing();
			Vector<B> roundedInteger = extension./* ringOfIntegers(). */asVector(this);
			List<Ext<B>> integralBasis = extension/* .ringOfIntegers() */.getModuleGenerators();
			StringBuilder build = new StringBuilder();
			boolean first = true;
			for (int i = 0; i < integralBasis.size(); i++) {
				B integer = lowAccuracyBase.round(roundedInteger.get(i + 1), lowAccuracyBase.getAccuracy());
				if (integer.equals(lowAccuracyBase.zero())) {
					continue;
				}
				if (first) {
					first = false;
				} else {
					build.append(" + ");
				}
				String coefficient;
				if (!integer.equals(lowAccuracyBase.one())) {
					coefficient = integer.toString();
				} else {
					coefficient = "";
				}
				String basis;
				if (!integralBasis.get(i).equals(extension.one())) {
					UnivariatePolynomial<B> roundedPolynomial = univariatePolynomialRing
							.getEmbedding(integralBasis.get(i).asPolynomial, new MathMap<>() {
								@Override
								public B evaluate(B t) {
									return lowAccuracyBase.round(t, lowAccuracyBase.getAccuracy());
								}
							});
					basis = roundedPolynomial.toString(/* "α", */true);
				} else {
					basis = "";
				}
				if (coefficient.equals("") && basis.equals("")) {
					build.append("1");
				} else if (coefficient.equals("")) {
					build.append(basis);
				} else if (basis.equals("")) {
					build.append(coefficient);
				} else {
					if (coefficient.contains(" ")) {
						build.append("(" + coefficient + ")");
					} else {
						build.append(coefficient);
					}
					build.append("*");
					if (basis.contains(" ")) {
						build.append("(" + basis + ")");
					} else {
						build.append(basis);
					}
				}
			}
			return build.toString();
		}

	}

	private DiscreteValuationField<B, S> baseField;
	private UnivariatePolynomial<B> minimalPolynomial;
	private Extension<S, R, RE, RFE> trivialReductionExtension;
	private OkutsuType<B, S, R, RE, RFE> type;
	private RFE residueField;
	private Ext<B> uniformizer;
	private List<Ext<B>> integralBasis;
	private Matrix<B> fromIntegralBasisBaseChange;
	private Matrix<B> toIntegralBasisBaseChange;
	private LocalRingExtension<B, S, Ext<B>, CompleteDVRExtension<B, S, R, RE, RFE>, R, RE, RFE> ringOfIntegers;
	private List<Value> integralBasisValues;
	private OtherVersion<Ext<B>, ?, RE, ?> exact;
	private int accuracy;

	public static <B extends Element<B>, S extends Element<S>> CompleteDVRExtension<B, S, ?, ?, ?> getCompleteDVRExtension(
			DiscreteValuationField<B, S> baseField) {
		return getCompleteDVRExtension(baseField,
				baseField.residueField().getExtension(baseField.residueField().getUnivariatePolynomialRing().getVar()));
	}

	public static <B extends Element<B>, S extends Element<S>, R extends Element<R>, RE extends AlgebraicExtensionElement<R, RE>, RFE extends FieldExtension<R, RE, RFE>> CompleteDVRExtension<B, S, R, RE, RFE> getCompleteDVRExtension(
			DiscreteValuationField<B, S> baseField, Extension<S, R, RE, RFE> trivialReductionExtension) {
		return new CompleteDVRExtension<B, S, R, RE, RFE>(baseField, trivialReductionExtension,
				baseField.getAccuracy());
	}

	public static <B extends Element<B>, S extends Element<S>> CompleteDVRExtension<B, S, ?, ?, ?> getCompleteDVRExtension(
			UnivariatePolynomial<B> minimalPolynomial, DiscreteValuationField<B, S> baseField) {
		return getCompleteDVRExtension(minimalPolynomial, baseField,
				baseField.residueField().getExtension(baseField.residueField().getUnivariatePolynomialRing().getVar()));
	}

	public static <B extends Element<B>, S extends Element<S>, R extends Element<R>, RE extends AlgebraicExtensionElement<R, RE>, RFE extends FieldExtension<R, RE, RFE>> CompleteDVRExtension<B, S, R, RE, RFE> getCompleteDVRExtension(
			UnivariatePolynomial<B> minimalPolynomial, DiscreteValuationField<B, S> baseField,
			Extension<S, R, RE, RFE> trivialReductionExtension) {
		TheMontesResult<B, S, R, RE, RFE> theMontes = baseField.ringOfIntegers().theMontesAlgorithm(minimalPolynomial,
				trivialReductionExtension);
		if (theMontes.getTypes().size() != 1) {
			throw new ArithmeticException("Polynomial not irreducible (at this accuracy)!");
		}
		OkutsuType<B, S, R, RE, RFE> type = theMontes.getTypes().get(0).updateLeaf(minimalPolynomial);
		int accuracy = baseField.getAccuracy();
		int initialAccuracy = type.previousLevel().precision().value();
		int requiredDigitAccuracy;
		if (type.quality().isInfinite()) {
			requiredDigitAccuracy = MiscAlgorithms.DivRoundUp(
					2 * type.ramificationIndex() * (Math.max(0, accuracy - initialAccuracy) + 1),
					type.ramificationIndex());
		} else {
			requiredDigitAccuracy = MiscAlgorithms
					.DivRoundUp(2 * type.ramificationIndex() * (Math.max(0, accuracy - initialAccuracy) + 1)
							+ type.quality().value(), type.ramificationIndex())
					+ 4 * type.exponent();
		}
		DiscreteValuationField<B, S> highAccuracyBaseField = baseField.withAccuracy(requiredDigitAccuracy);
		UnivariatePolynomialRing<B> highAccuracyPolynomialRing = highAccuracyBaseField.getUnivariatePolynomialRing();
		UnivariatePolynomial<B> highAccuracyMinimalPolynomial = highAccuracyPolynomialRing
				.getEmbedding(minimalPolynomial, new MathMap<>() {
					@Override
					public B evaluate(B t) {
						return highAccuracyBaseField.round(t, highAccuracyBaseField.getAccuracy() + 1);
					}
				});
		return new CompleteDVRExtension<>(highAccuracyMinimalPolynomial, highAccuracyBaseField,
				trivialReductionExtension, accuracy);
	}

	public static <E extends Element<E>, EE extends AlgebraicExtensionElement<E, EE>, EFE extends FieldExtension<E, EE, EFE> & DiscreteValuationField<EE, RE>, S extends Element<S>, R extends Element<R>, RE extends AlgebraicExtensionElement<R, RE>, RFE extends FieldExtension<R, RE, RFE>> CompleteDVRExtension<?, S, R, RE, RFE> getCompleteDVRExtension(
			EFE exact, OkutsuType<E, S, R, RE, RFE> type, int accuracy,
			Extension<S, R, RE, RFE> trivialReductionExtension) {
		int initialAccuracy = type.previousLevel().precision().value();
		int requiredDigitAccuracy;
		if (type.quality().isInfinite()) {
			requiredDigitAccuracy = MiscAlgorithms.DivRoundUp(
					2 * type.ramificationIndex() * Math.max(0, accuracy - initialAccuracy), type.ramificationIndex());
		} else {
			requiredDigitAccuracy = MiscAlgorithms.DivRoundUp(
					2 * type.ramificationIndex() * Math.max(0, accuracy - initialAccuracy) + type.quality().value(),
					type.ramificationIndex())
			/* + 4 * type.exponent() */;
		}
		@SuppressWarnings("unchecked")
		DiscreteValuationField<E, S> exactBaseField = (DiscreteValuationField<E, S>) exact.getBaseField();
		type = exactBaseField.ringOfIntegers().singleFactorLifting(type, accuracy);
		return getCompleteDVRExtension(exact, type, accuracy, exactBaseField.complete(requiredDigitAccuracy),
				trivialReductionExtension);
	}

	private static <E extends Element<E>, EE extends AlgebraicExtensionElement<E, EE>, EFE extends FieldExtension<E, EE, EFE> & DiscreteValuationField<EE, RE>, B extends Element<B>, S extends Element<S>, R extends Element<R>, RE extends AlgebraicExtensionElement<R, RE>, RFE extends FieldExtension<R, RE, RFE>, DVRB extends DiscreteValuationField<B, S>> CompleteDVRExtension<B, S, R, RE, RFE> getCompleteDVRExtension(
			EFE exact, OkutsuType<E, S, R, RE, RFE> type, int accuracy, OtherVersion<E, B, S, DVRB> completedBaseField,
			Extension<S, R, RE, RFE> trivialReductionExtension) {
		DiscreteValuationField<B, S> baseField = completedBaseField.getField();
		UnivariatePolynomialRing<B> polynomialRing = baseField.getUnivariatePolynomialRing();
		UnivariatePolynomial<B> minimalPolynomial = polynomialRing.getEmbedding(type.representative(),
				completedBaseField.getRetraction());
		return new CompleteDVRExtension<>(minimalPolynomial, baseField, trivialReductionExtension, accuracy);
	}

	private CompleteDVRExtension(DiscreteValuationField<B, S> baseField,
			Extension<S, R, RE, RFE> trivialReductionExtension, int accuracy) {
		super(baseField);
		this.accuracy = accuracy;
		this.baseField = baseField;
		this.minimalPolynomial = baseField.getUnivariatePolynomialRing().getVar();
		init(trivialReductionExtension);
	}

	private CompleteDVRExtension(UnivariatePolynomial<B> minimalPolynomial, DiscreteValuationField<B, S> baseField,
			Extension<S, R, RE, RFE> trivialReductionExtension, int accuracy) {
		super(minimalPolynomial, baseField, "x");
		this.accuracy = accuracy;
		this.baseField = baseField;
		this.minimalPolynomial = minimalPolynomial;
		init(trivialReductionExtension);
	}

	private void init(Extension<S, R, RE, RFE> trivialReductionExtension) {
		TheMontesResult<B, S, R, RE, RFE> theMontes = baseField.ringOfIntegers().theMontesAlgorithm(minimalPolynomial,
				trivialReductionExtension);
		if (theMontes.getTypes().size() != 1) {
			throw new ArithmeticException("Polynomial not irreducible (at this accuracy)!");
		}
		this.trivialReductionExtension = trivialReductionExtension;
		this.type = theMontes.getTypes().get(0).updateLeaf(minimalPolynomial);
		this.residueField = type.reduction().extension();
		this.uniformizer = fromSmallDegreePolynomial(type.lift(residueField().one(), 1));
		List<UnivariatePolynomial<B>> polynomialBasis = baseField.ringOfIntegers().triagonalizeIntegralBasis(
				minimalPolynomial, baseField.ringOfIntegers().integralBasis(minimalPolynomial, theMontes, true));
		this.integralBasis = new ArrayList<>();
		List<Vector<B>> basisVectors = new ArrayList<>();
		integralBasisValues = new ArrayList<>();
		for (UnivariatePolynomial<B> polynomial : polynomialBasis) {
			Ext<B> vector = fromSmallDegreePolynomial(polynomial);
			integralBasisValues.add(valuation(vector));
			integralBasis.add(vector);
			basisVectors.add(asVector(vector));
		}
		this.fromIntegralBasisBaseChange = Matrix.fromColumns(basisVectors);
		this.toIntegralBasisBaseChange = matrixAlgebra().invertUpperTriangleMatrix(fromIntegralBasisBaseChange);
		this.ringOfIntegers = new LocalRingExtension<>(this, baseField, residueField, integralBasis,
				toIntegralBasisBaseChange);
	}

	public int ramificationIndex() {
		return type.ramificationIndex();
	}

	public int residueDegree() {
		return type.residueDegree();
	}

	@Override
	public Real value(Ext<B> t) {
		return Reals.r(1024).positiveRoot(baseField.value(norm(t)), degree());
	}

	@Override
	public Reals getReals() {
		return Reals.r(1024);
	}

	@Override
	public boolean isComplete() {
		return true;
	}

	@Override
	public Ext<B> inverse(Ext<B> t, int accuracy) {
		if (accuracy > getAccuracy()) {
			return withAccuracy(accuracy).inverse(t, accuracy);
		}
		Value value = valuation(t);
		if (value.equals(Value.INFINITY)) {
			throw new ArithmeticException("Division by 0!");
		}
		if (value.equals(Value.ZERO)) {
			return fromSmallDegreePolynomial(baseField.ringOfIntegers().invertInType(t.asPolynomial(), type, accuracy));
		}
		if (value.compareTo(Value.ZERO) < 0) {
			Ext<B> uniformizerPower = fromSmallDegreePolynomial(type.lift(residueField().one(), -value.value()));
			return multiply(inverse(multiply(t, uniformizerPower), accuracy), uniformizerPower);
		}
		Ext<B> uniformizerPower = fromSmallDegreePolynomial(type.lift(residueField().one(), -value.value()));
		return multiply(inverse(multiply(t, uniformizerPower), accuracy + value.value()), uniformizerPower);
	}

	@Override
	public Ext<B> inverse(Ext<B> t) {
		return inverse(t, getAccuracy());
	}

	@Override
	public Ext<B> divide(Ext<B> dividend, Ext<B> divisor) {
		if (dividend.equals(zero())) {
			return zero();
		}
		Value dividendValue = valuation(dividend);
		Value divisorValue = valuation(divisor);
		int accuracy = getAccuracy() - Math.min(0, dividendValue.value()) - Math.min(0, divisorValue.value());
		return multiply(dividend, inverse(divisor, accuracy));
	}

	@Override
	public Ext<B> negative(Ext<B> t, int accuracy) {
		Vector<B> asVector = ringOfIntegers().asVector(t);
		List<B> negative = new ArrayList<>();
		for (B b : asVector.asList()) {
			negative.add(baseField.negative(b, accuracy));
		}
		return ringOfIntegers().fromVector(new Vector<>(negative));
	}

	@Override
	public Value valuation(Ext<B> t) {
		return type.valuation(t.asPolynomial());
	}

	@Override
	public Ext<B> uniformizer() {
		return uniformizer;
	}

	@Override
	public Field<RE> residueField() {
		return residueField;
	}

	@Override
	public RE reduceInteger(Ext<B> t) {
		Value v = valuation(t);
		if (v.compareTo(Value.ZERO) > 0) {
			return residueField().zero();
		} else if (v.compareTo(Value.ZERO) < 0) {
			throw new ArithmeticException("not an integer!");
		}
		return type.reduce(t.asPolynomial());
	}

	@Override
	public Ext<B> upToUniformizerPower(Ext<B> t) {
		return fromSmallDegreePolynomial(type.lift(type.reduce(t.asPolynomial()), 0));
	}

	@Override
	public Ext<B> liftToInteger(RE s) {
		return fromSmallDegreePolynomial(type.lift(s, 0));
	}

	@Override
	public Ext<B> round(Ext<B> t, int accuracy) {
		Vector<B> asIntegerVector = ringOfIntegers().asVector(t);
		int ramification = type.ramificationIndex();
		List<B> roundedVector = new ArrayList<>();
		for (int i = 0; i < asIntegerVector.dimension(); i++) {
			B coefficient = asIntegerVector.get(i + 1);
			roundedVector.add(baseField.round(coefficient,
					MiscAlgorithms.DivRoundUp(accuracy - integralBasisValues.get(i).value(), ramification)));
		}
		return ringOfIntegers().fromVector(new Vector<>(roundedVector));
	}

	@Override
	public int getAccuracy() {
		return accuracy;
	}

	@Override
	public DiscreteValuationField<B, S> getBaseField() {
		return baseField;
	}

	@Override
	public CompleteDVRExtension<B, S, R, RE, RFE> withAccuracy(int accuracy) {
		return getCompleteDVRExtension(minimalPolynomial, baseField.withAccuracy(accuracy), trivialReductionExtension);
	}

	public void setExact(OtherVersion<Ext<B>, ?, RE, ?> exact) {
		if (this.exact != null) {
			return;
		}
		this.exact = exact;
	}

	@Override
	public OtherVersion<Ext<B>, ?, RE, ?> exact() {
		if (exact == null) {
			exact = exact(baseField.exact());
		}
		return exact;
	}

	@SuppressWarnings("unchecked")
	private <T extends Element<T>, Ex extends DiscreteValuationField<T, S>> OtherVersion<Ext<B>, ?, RE, ?> exact(
			OtherVersion<B, T, S, Ex> exactBase) {
		UnivariatePolynomialRing<T> exactBasePolynomialRing = exactBase.getField().getUnivariatePolynomialRing();
		return (OtherVersion<Ext<B>, ?, RE, ?>) exact(exactBase, exactBase.getField().getUniqueExtension(
				exactBasePolynomialRing.getEmbedding(minimalPolynomial, exactBase.getRetraction())));
	}

	
	private <T extends Element<T>,SR extends Element<SR>, Ex extends DiscreteValuationField<T, SR>,SE extends Element<SE>, E extends Element<E>, EE extends AlgebraicExtensionElement<E, EE>,REx extends Element<REx>, REEx extends AlgebraicExtensionElement<REx, REEx>, RFEx extends FieldExtension<REx, REEx, RFEx>, EFE extends DiscreteValuationFieldExtension<E, SE, EE, EFE, REx, REEx, RFEx>> OtherVersion<Ext<B>, EE, REEx, EFE> exact(
			OtherVersion<B, T, SR, Ex> exactBase, ExtensionOfDiscreteValuationField<T, SR, E, SE, EE, EFE, REx, REEx, RFEx> extension) {
		UnivariatePolynomialRing<T> exactBasePolynomialRing = exactBase.getField().getUnivariatePolynomialRing();
		EFE ext = extension.getExtension();
		UnivariatePolynomialRing<EE> exactExtensionPolynomialRing = ext.getUnivariatePolynomialRing();
		return new OtherVersion<Ext<B>, EE, REEx, EFE>(ext, new MathMap<>() {

			@Override
			public EE evaluate(Ext<B> t) {
				return exactExtensionPolynomialRing
						.evaluate(
								exactExtensionPolynomialRing.getEmbedding(exactBasePolynomialRing.getEmbedding(
										t.asPolynomial(), exactBase.getRetraction()), extension.getEmbeddingMap()),
								Collections.singletonList(ext.alpha()));
			}
		}, new MathMap<>() {

			@Override
			public Ext<B> evaluate(EE t) {
				Vector<T> asVector = extension.getAsVectorMap().evaluate(t);
				List<B> asRoundedVector = new ArrayList<>();
				for (T c : asVector.asList()) {
					asRoundedVector.add(exactBase.getEmbedding().evaluate(c));
				}
				return fromVector(new Vector<>(asRoundedVector));
			}
		});
	}

	@Override
	public OtherVersion<Ext<B>, Ext<B>, RE, CompleteDVRExtension<B, S, R, RE, RFE>> complete(int accuracy) {
		return new OtherVersion<>(withAccuracy(accuracy), new Identity<>(), new Identity<>());
	}

	@Override
	public Ext<B> getRandomInteger() {
		List<B> asVector = new ArrayList<>();
		for (int i = 0; i < degree(); i++) {
			asVector.add(baseField.getRandomInteger());
		}
		return ringOfIntegers().fromVector(new Vector<>(asVector));
	}

	@Override
	public LocalRingExtension<B, S, Ext<B>, CompleteDVRExtension<B, S, R, RE, RFE>, R, RE, RFE> ringOfIntegers() {
		return ringOfIntegers;
	}

	@Override
	protected Ext<B> fromSmallDegreePolynomial(UnivariatePolynomial<B> polynomial) {
		return new Ext<>(this, polynomial);
	}

	@Override
	public CompleteDVRExtension<B, S, R, RE, RFE> makeExtension(UnivariatePolynomial<B> minimalPolynomial) {
		return getCompleteDVRExtension(minimalPolynomial, baseField, trivialReductionExtension);
	}

	@Override
	public ExtensionOfDiscreteValuationField<Ext<B>, RE, B, S, Ext<B>, CompleteDVRExtension<B, S, R, RE, RFE>, R, RE, RFE> getUniqueExtension(
			UnivariatePolynomial<Ext<B>> minimalPolynomial) {
		Extension<Ext<B>, B, Ext<B>, CompleteDVRExtension<B, S, R, RE, RFE>> extension = getExtension(
				minimalPolynomial);
		return new ExtensionOfDiscreteValuationField<>(this, extension.extension(), extension.embeddingMap(),
				extension.asVectorMap());
	}

	@Override
	protected CompleteDVRExtension<B, S, R, RE, RFE> asExtensionType() {
		return this;
	}

	@Override
	public FactorizationResult<Polynomial<Ext<B>>, Ext<B>> factorization(UnivariatePolynomial<Ext<B>> t) {
		return ringOfIntegers().factorization(t, true, getAccuracy());
	}

}
