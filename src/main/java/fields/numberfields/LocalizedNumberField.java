package fields.numberfields;

import java.util.ArrayList;
import java.util.List;

import fields.exceptions.InfinityException;
import fields.finitefields.FiniteField;
import fields.finitefields.FiniteField.FFE;
import fields.finitefields.PrimeField.PFE;
import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.helper.AbstractFieldExtension;
import fields.integers.Integers.IntE;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.DiscreteValuationField;
import fields.interfaces.DiscreteValuationRing;
import fields.interfaces.DiscreteValuationRing.OkutsuType;
import fields.interfaces.MathMap;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.local.CompleteDVRExtension;
import fields.local.CompleteDVRExtension.Ext;
import fields.local.LocalRingExtension;
import fields.local.PAdicField;
import fields.local.PAdicField.PAdicNumber;
import fields.local.Value;
import fields.numberfields.NumberField.NFE;
import fields.numberfields.NumberFieldIntegers.NumberFieldIdeal;
import fields.vectors.pivot.PivotStrategy;
import fields.vectors.pivot.ValuationPivotStrategy;
import util.Identity;

public class LocalizedNumberField extends AbstractFieldExtension<Fraction, NFE, LocalizedNumberField>
		implements DiscreteValuationField<NFE, FFE> {
	private NumberField field;
	private IntE prime;
	private OkutsuType<Fraction, PFE, PFE, FFE, FiniteField> type;
	private NumberFieldIdeal ideal;
	private LocalRingExtension<Fraction, PFE, NFE, LocalizedNumberField, PFE, FFE, FiniteField> localRing;
	private FiniteField reduction;
	private NFE uniformizer;
	private Reals r = Reals.r(1024);

	LocalizedNumberField(NumberField field, NumberFieldIdeal ideal,
			DiscreteValuationRing<Fraction, PFE> localizedIntegers) {
		super(field.minimalPolynomial(), localizedIntegers.localField());
		this.field = field;
		this.ideal = ideal;
		this.prime = ideal.prime();
		this.type = ideal.type();
		this.reduction = type.reduction().extension();
		this.uniformizer = field.fromPolynomial(type.lift(reduction.one(), 1));
		this.localRing = new LocalRingExtension<>(this, localizedIntegers.localField(), type,
				field.maximalOrder().getModuleGenerators(), field.maximalOrder().toIntegralBasisBaseChange());
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

	@Override
	public NFE uniformizer() {
		return uniformizer;
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
		return field.fromPolynomial(type.lift(s, 0));
	}

	@Override
	public NFE round(NFE t, int accuracy) {
		return ideal.round(t, accuracy);
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
	public OtherVersion<NFE, Ext<PAdicNumber>, FFE, CompleteDVRExtension<PAdicNumber, PFE, PFE, FFE, FiniteField>> complete(
			int accuracy) {
		UnivariatePolynomialRing<NFE> nfPolynomials = field.getUnivariatePolynomialRing();
		PAdicField base = new PAdicField(prime.getValue(), accuracy);
		UnivariatePolynomialRing<PAdicNumber> polynomials = base.getUnivariatePolynomialRing();
		Extension<PAdicNumber, PAdicNumber, Ext<PAdicNumber>, CompleteDVRExtension<PAdicNumber, PFE, PFE, FFE, FiniteField>> extension = base
				.getExtension(polynomials.getEmbedding(type.representative(), base.fromRationalMap()));
		UnivariatePolynomialRing<Ext<PAdicNumber>> extensionPolynomials = extension.extension()
				.getUnivariatePolynomialRing();
		Ext<PAdicNumber> alpha = extension.extension().alpha();
		return new OtherVersion<>(extension.extension(), new MathMap<NFE, Ext<PAdicNumber>>() {

			@SuppressWarnings("unchecked")
			@Override
			public Ext<PAdicNumber> evaluate(NFE t) {
				return extensionPolynomials.evaluate(extensionPolynomials.getEmbedding(
						polynomials.getEmbedding(t.asPolynomial(), base.fromRationalMap()), extension.embeddingMap()),
						alpha);
			}
		}, new MathMap<Ext<PAdicNumber>, NFE>() {

			@Override
			public NFE evaluate(Ext<PAdicNumber> t) {
				return nfPolynomials.evaluate(
						nfPolynomials.getEmbedding(Rationals.q().getUnivariatePolynomialRing()
								.getEmbedding(t.asPolynomial(), base.toRationalMap()), field.getEmbeddingMap()),
						field.alpha());
			}
		});
	}

	@Override
	public DiscreteValuationFieldExtension<NFE, FFE, Fraction, NFE, FFE, LocalizedNumberField> getUniqueExtension(
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
		return new DiscreteValuationFieldExtension<>(this, maximalOrder.localizeAndQuotient(newIdeal),
				extension.embeddingMap(), extension.asVectorMap());
	}

}
