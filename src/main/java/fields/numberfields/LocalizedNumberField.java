package fields.numberfields;

import java.util.ArrayList;
import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;

import fields.exceptions.InfinityException;
import fields.finitefields.FiniteField;
import fields.finitefields.FiniteField.FFE;
import fields.finitefields.PrimeField.PFE;
import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.helper.AbstractFieldExtension;
import fields.integers.Integers.IntE;
import fields.integers.LocalizedFractions;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.DiscreteValuationField;
import fields.interfaces.DiscreteValuationFieldExtension;
import fields.interfaces.DiscreteValuationRing;
import fields.interfaces.DiscreteValuationRing.OkutsuType;
import fields.interfaces.UnivariatePolynomial;
import fields.local.LocalRingExtension;
import fields.local.Value;
import fields.numberfields.CompletedNumberField.Ext;
import fields.numberfields.NumberField.NFE;
import fields.numberfields.NumberFieldIntegers.NumberFieldIdeal;
import fields.vectors.pivot.PivotStrategy;
import fields.vectors.pivot.ValuationPivotStrategy;
import util.Identity;
import util.MiscAlgorithms;

public class LocalizedNumberField extends AbstractFieldExtension<Fraction, NFE, LocalizedNumberField>
		implements DiscreteValuationFieldExtension<Fraction, PFE, NFE, LocalizedNumberField, PFE, FFE, FiniteField> {
	private NumberField field;
	private IntE prime;
	private OkutsuType<Fraction, PFE, PFE, FFE, FiniteField> type;
	private NumberFieldIdeal ideal;
	private LocalRingExtension<Fraction, PFE, NFE, LocalizedNumberField, PFE, FFE, FiniteField> localRing;
	private FiniteField reduction;
	private NFE uniformizer;
	private Reals r;
	private SortedMap<Integer, OtherVersion<NFE, Ext, FFE, CompletedNumberField>> complete;

	LocalizedNumberField(NumberField field, NumberFieldIdeal ideal,
			DiscreteValuationRing<Fraction, PFE> localizedIntegers) {
		super(field.minimalPolynomial(), localizedIntegers.localField(), field.getVariableName());
		this.r = field.minkowskiEmbeddingSpace().getValueField().getReals();
		this.field = field;
		this.ideal = ideal;
		this.prime = ideal.prime();
		this.type = ideal.type();
		this.reduction = type.reduction().extension();
		this.uniformizer =/* ideal.isPrincipal() ? ideal.principalGenerator() :*/ ideal.uniformizer();// field.fromPolynomial(type.lift(reduction.one(), 1));
		this.localRing = new LocalRingExtension<>(this, localizedIntegers.localField(), reduction,
				field.maximalOrder().getModuleGenerators(), field.maximalOrder().toIntegralBasisBaseChange());
		this.complete = new TreeMap<>();
	}

	@Override
	public Real value(NFE t) {
		if (t.equals(zero())) {
			return r.zero();
		}
		return r.exp(
				r.divide(r.multiply(r.log(r.getInteger(prime)), r.getInteger(type.valuation(t.asPolynomial()).value())),
						r.getInteger(type.ramificationIndex())));
	}

	@Override
	public Reals getReals() {
		return r;
	}
	
	@Override
	public LocalizedFractions getBaseField() {
		return Rationals.q().withValuation(prime);
	}

	@Override
	public NFE inverse(NFE t, int accuracy) {
		return inverse(t);
	}

	@Override
	public NFE negative(NFE t, int accuracy) {
		return negative(t);
	}

	@Override
	public PivotStrategy<NFE> preferredPivotStrategy() {
		return new ValuationPivotStrategy<>(this.valuation());
	}

	@Override
	public boolean isInteger(NFE t) {
		return valuation(t).compareTo(Value.ZERO) >= 0;
	}

	@Override
	public Value valuation(NFE t) {
		return type.valuation(t.asPolynomial());
	}

	public int ramificationIndex() {
		return type().ramificationIndex();
	}

	public int residueDegree() {
		return type.residueDegree();
	}

	public int completedDegree() {
		return ramificationIndex() * residueDegree();
	}

	@Override
	public NFE uniformizer() {
		return uniformizer;
	}

	OkutsuType<Fraction, PFE, PFE, FFE, FiniteField> type() {
		return type;
	}

	NumberField getNumberField() {
		return field;
	}

	NumberFieldIdeal ideal() {
		return ideal;
	}

	@Override
	public FiniteField residueField() {
		return reduction;
	}

	@Override
	public FFE reduceInteger(NFE t) {
		if (!type.valuation(t.asPolynomial()).equals(Value.ZERO)) {
			return reduction.zero();
		}
		return type.reduce(t.asPolynomial());
	}

	@Override
	public NFE upToUniformizerPower(NFE t) {
		if (t.equals(zero())) {
			return zero();
		}
		return liftToInteger(type.reduce(t.asPolynomial()));
	}

	@Override
	public NFE liftToInteger(FFE s) {
		return round(field.fromPolynomial(type.lift(s, 0)), 1);
	}

	@Override
	public NFE round(NFE t, int accuracy) {
		OtherVersion<NFE, Ext, FFE, CompletedNumberField> complete;
		if (!this.complete.isEmpty() && this.complete.lastKey() > accuracy) {
			complete = complete(this.complete.lastKey());
		} else {
			complete = complete(accuracy + 1);
		}
		return complete.getEmbedding()
				.evaluate(complete.getField().round(complete.getRetraction().evaluate(t), accuracy));
	}

	@Override
	public int getAccuracy() {
		throw new InfinityException();
	}

	@Override
	public DiscreteValuationField<NFE, FFE> withAccuracy(int accuracy) {
		return this;
	}

	@Override
	public NFE getRandomInteger() {
		return field.maximalOrder().getRandomElement();
	}

	@Override
	public LocalRingExtension<Fraction, PFE, NFE, LocalizedNumberField, PFE, FFE, FiniteField> ringOfIntegers() {
		return localRing;
	}

	@Override
	protected NFE fromSmallDegreePolynomial(UnivariatePolynomial<Fraction> polynomial) {
		if (field == null) {
			return new NFE(null, polynomial);
		}
		return field.fromSmallDegreePolynomial(polynomial);
	}

	@Override
	public LocalizedNumberField makeExtension(UnivariatePolynomial<Fraction> minimalPolynomial) {
		throw new UnsupportedOperationException("LocalizedNumberField requires additional information!");
	}

	@Override
	protected LocalizedNumberField asExtensionType() {
		return this;
	}

	@Override
	public boolean isComplete() {
		return false;
	}

	@Override
	public OtherVersion<NFE, NFE, FFE, LocalizedNumberField> exact() {
		return new OtherVersion<>(this, new Identity<>(), new Identity<>());
	}

	@Override
	public OtherVersion<NFE, Ext, FFE, CompletedNumberField> complete(int accuracy) {
		if (accuracy < 4) {
			accuracy = 4;
		}
		int notRounded = accuracy;
		accuracy = 1;
		while (notRounded != 1) {
			notRounded = MiscAlgorithms.DivRoundUp(notRounded, 2);
			accuracy *= 2;
		}
		if (!complete.containsKey(accuracy)) {
//			PrimeField reductionBase = PrimeField.getPrimeField(prime);
//			@SuppressWarnings("unchecked")
//			CompleteDVRExtension<PAdicNumber, PFE, PFE, FFE, FiniteField> extension = (CompleteDVRExtension<PAdicNumber, PFE, PFE, FFE, FiniteField>) CompleteDVRExtension
//					.getCompleteDVRExtension(this, type, accuracy,
//							reductionBase.getExtension(reductionBase.getUnivariatePolynomialRing().getVar()));
//			PAdicField base = (PAdicField) extension.getBaseField();
//			UnivariatePolynomialRing<PAdicNumber> polynomials = base.getUnivariatePolynomialRing();
//			UnivariatePolynomialRing<Ext<PAdicNumber>> extensionPolynomials = extension.getUnivariatePolynomialRing();
//			Ext<PAdicNumber> alpha = extension.alpha();
//			MathMap<NFE, Ext<PAdicNumber>> embeddingMap = new MathMap<NFE, Ext<PAdicNumber>>() {
//
//				@SuppressWarnings("unchecked")
//				@Override
//				public Ext<PAdicNumber> evaluate(NFE t) {
//					return extensionPolynomials.evaluate(extensionPolynomials.getEmbedding(
//							polynomials.getEmbedding(t.asPolynomial(), base.fromRationalMap()),
//							extension.getEmbeddingMap()), alpha);
//				}
//			};
//			List<Vector<IntE>> asIntegerVectors = new ArrayList<>();
//			Integers z = Integers.z();
//			for (NFE integral : field.maximalOrder().getModuleGenerators()) {
//				Vector<PAdicNumber> asVector = extension.asVector(embeddingMap.evaluate(integral));
//				asIntegerVectors.add(Vector.mapVector(base.roundToIntegerMap(base.getAccuracy()), asVector));
//			}
//			SmallestIntegerSolutionPreparation preparation = z.prepareSmallestIntegerSolution(asIntegerVectors,
//					z.power(prime,MiscAlgorithms.DivRoundDown(accuracy, type.ramificationIndex())), 1.0);
//			MathMap<Ext<PAdicNumber>, NFE> retractionMap = new MathMap<>() {
//				@Override
//				public NFE evaluate(Ext<PAdicNumber> t) {
//					//t = extension.round(t, finalAccuracy);
//					Vector<PAdicNumber> asVector = extension.asVector(t);
//					Value valuation = Value.INFINITY;
//					for (PAdicNumber coeff : asVector.asList()) {
//						valuation = valuation.min(base.valuation(coeff));
//					}
//					PAdicField accurateBase = base;
//					if (valuation.compareTo(Value.ZERO) < 0) {
//						accurateBase = accurateBase.withAccuracy(accurateBase.getAccuracy() - valuation.value());
//						List<PAdicNumber> integers = new ArrayList<>();
//						PAdicNumber uniformizerPower = accurateBase.power(accurateBase.uniformizer(),
//								-valuation.value());
//						for (PAdicNumber coeff : asVector.asList()) {
//							integers.add(accurateBase.multiply(coeff, uniformizerPower));
//						}
//						asVector = new Vector<>(integers);
//					}
//					Vector<IntE> asIntegerVector = Vector.mapVector(accurateBase.roundToIntegerMap(accurateBase.getAccuracy()),
//							asVector);
//					Vector<IntE> inIntegralBasis = z.smallestIntegerSolution(asIntegerVector, preparation);
//					NFE result = field.maximalOrder().fromVector(inIntegralBasis);
//					if (valuation.compareTo(Value.ZERO) < 0) {
//						NFE uniformizerPowerInverse = field.power(field.getEmbedding(prime), valuation.value());
//						result = field.multiply(uniformizerPowerInverse, result);
//					}
//					return result;
//				}
//			};
//			extension.setExact(new OtherVersion<>(this, retractionMap, embeddingMap));
			CompletedNumberField completed = CompletedNumberField.getCompletedNumberField(this, accuracy);
			complete.put(accuracy,
					new OtherVersion<>(completed, completed.exact().getEmbedding(), completed.exact().getRetraction()));
		}
		return complete.get(accuracy);
	}

	@Override
	public ExtensionOfDiscreteValuationField<NFE, FFE, Fraction, PFE, NFE, LocalizedNumberField,PFE, FFE, FiniteField> getUniqueExtension(
			UnivariatePolynomial<NFE> minimalPolynomial) {
		Extension<NFE, Fraction, NFE, NumberField> extension = field.getExtension(minimalPolynomial);
		NumberFieldIntegers maximalOrder = extension.extension().maximalOrder();
		List<NFE> generators = new ArrayList<>();
		for (NFE generator : ideal.generators()) {
			generators.add(extension.embeddingMap().evaluate(generator));
		}
		NumberFieldIdeal newIdeal = maximalOrder.radical(maximalOrder.getIdeal(generators));
		if (!newIdeal.isMaximal()) {
			throw new ArithmeticException("Field extension is split!");
		}
		return new ExtensionOfDiscreteValuationField<>(this, maximalOrder.localizeAndQuotient(newIdeal),
				extension.embeddingMap(), extension.asVectorMap());
	}

}
