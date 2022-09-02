package fields.local;

import java.io.IOException;
import java.math.BigInteger;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

import fields.exceptions.InfinityException;
import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.helper.AbstractField;
import fields.helper.TranscendentalFieldExtension;
import fields.helper.TranscendentalFieldExtension.TExt;
import fields.interfaces.DiscreteValuationField;
import fields.interfaces.DiscreteValuationRing;
import fields.interfaces.Element;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.UnivariatePolynomial;
import fields.polynomials.Monomial;
import fields.vectors.pivot.PivotStrategy;
import fields.vectors.pivot.ValuationPivotStrategy;
import util.PeekableReader;
import varieties.SimpleFunctionField;
import varieties.SimpleFunctionField.SimpleRationalFunction;

public class TranscendentalDVRExtension<T extends Element<T>, S extends Element<S>> extends AbstractField<TExt<T>>
		implements DiscreteValuationField<TExt<T>, TExt<S>> {
	private TranscendentalFieldExtension<T> extension;
	private DiscreteValuationField<T, S> base;
	private LocalRingImplementation<TExt<T>, TExt<S>> localRing;
	private TranscendentalFieldExtension<S> reduction;

	public TranscendentalDVRExtension(DiscreteValuationField<T, S> base, String[] variableNames) {
		this.extension = new TranscendentalFieldExtension<>(base, variableNames);
		this.base = base;
		this.reduction = new TranscendentalFieldExtension<>(base.residueField(),
				extension.polynomialRing().getVariableNames());
		this.localRing = new LocalRingImplementation<>(this, base.ringOfIntegers().toString() + "("
				+ String.join(", ", extension.polynomialRing().getVariableNames()) + ")");
	}

	public String toString() {
		return extension.toString();
	}

	public TExt<T> parse(PeekableReader reader) throws IOException {
		return extension.parse(reader);
	}

	public boolean equals(Object obj) {
		return extension.equals(obj);
	}

	public TExt<T> zero() {
		return extension.zero();
	}

	public TExt<T> add(TExt<T> s1, TExt<T> s2) {
		return extension.add(s1, s2);
	}

	public TExt<T> negative(TExt<T> s) {
		return extension.negative(s);
	}

	public Exactness exactness() {
		return extension.exactness();
	}

	public TExt<T> getRandomElement() {
		return extension.getRandomElement();
	}

	public boolean isFinite() {
		return extension.isFinite();
	}

	public BigInteger getNumberOfElements() throws InfinityException {
		return extension.getNumberOfElements();
	}

	public Iterator<TExt<T>> iterator() {
		return extension.iterator();
	}

	public TExt<T> one() {
		return extension.one();
	}

	public BigInteger characteristic() {
		return base.characteristic();
	}

	public TExt<T> multiply(TExt<T> t1, TExt<T> t2) {
		return extension.multiply(t1, t2);
	}

	public TExt<T> inverse(TExt<T> t) {
		return extension.inverse(t);
	}

	public TExt<T> getEmbedding(TExt<T> t) {
		return t;
	}

	public TExt<T> getEmbedding(T t) {
		return extension.getEmbedding(t);
	}

	public TExt<T> getEmbedding(Polynomial<T> t) {
		return extension.getEmbedding(t);
	}

	public TExt<T> getVar(int i) {
		return extension.getVar(i);
	}

	public TExt<T> getVarPower(int i, int power) {
		return extension.getVarPower(i, power);
	}

	public TranscendentalDVRExtension<T, S> asDedekindRing() {
		return this;
	}

	public DiscreteValuationField<T, S> getField() {
		return base;
	}

	public int transcendenceDegree() {
		return extension.transcendenceDegree();
	}

	public PolynomialRing<T> polynomialRing() {
		return extension.polynomialRing();
	}

	public Extension<TExt<T>, TExt<T>, SimpleRationalFunction<T>, SimpleFunctionField<T>> getExtension(
			UnivariatePolynomial<TExt<T>> minimalPolynomial) {
		return extension.getExtension(minimalPolynomial);
	}

	public FactorizationResult<Polynomial<TExt<T>>, TExt<T>> factorization(UnivariatePolynomial<TExt<T>> t) {
		return extension.factorization(t);
	}

	public PivotStrategy<TExt<T>> preferredPivotStrategy() {
		return new ValuationPivotStrategy<>(ringOfIntegers().valuation());
	}

	public TExt<T> characteristicRoot(TExt<T> t, int power) {
		return extension.characteristicRoot(t, power);
	}

	@Override
	public boolean hasCharacteristicRoot(TExt<T> t, int power) {
		return extension.hasCharacteristicRoot(t, power);
	}

	@Override
	public Real value(TExt<T> t) {
		Real result = getReals().zero();
		for (Monomial m : t.getNumerator().monomials()) {
			result = getReals().max(result, base.value(t.getNumerator().coefficient(m)));
		}
		return result;
	}

	@Override
	public Reals getReals() {
		return base.getReals();
	}

	@Override
	public boolean isComplete() {
		return base.isComplete();
	}

	@Override
	public TExt<T> inverse(TExt<T> t, int accuracy) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public TExt<T> negative(TExt<T> t, int accuracy) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Value valuation(TExt<T> t) {
		Value v = Value.INFINITY;
		for (Monomial m : t.getNumerator().monomials()) {
			v = v.min(base.valuation(t.getNumerator().coefficient(m)));
		}
		return v;
	}

	@Override
	public TExt<T> uniformizer() {
		return extension.getEmbedding(base.uniformizer());
	}

	@Override
	public TranscendentalFieldExtension<S> residueField() {
		return reduction;
	}

	@Override
	public TExt<S> reduceInteger(TExt<T> t) {
		Map<Monomial, S> numeratorMap = new TreeMap<>();
		for (Monomial m : t.getNumerator().monomials()) {
			numeratorMap.put(reduction.polynomialRing().getMonomial(m.exponents()),
					base.reduceInteger(t.getNumerator().coefficient(m)));
		}
		Map<Monomial, S> denominatorMap = new TreeMap<>();
		for (Monomial m : t.getDenominator().monomials()) {
			denominatorMap.put(reduction.polynomialRing().getMonomial(m.exponents()),
					base.reduceInteger(t.getDenominator().coefficient(m)));
		}
		return reduction.divide(reduction.getEmbedding(reduction.polynomialRing().getPolynomial(numeratorMap)),
				reduction.getEmbedding(reduction.polynomialRing().getPolynomial(denominatorMap)));
	}

	@Override
	public TExt<T> upToUniformizerPower(TExt<T> t) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public TExt<T> liftToInteger(TExt<S> s) {
		Map<Monomial, T> numeratorMap = new TreeMap<>();
		for (Monomial m : s.getNumerator().monomials()) {
			numeratorMap.put(extension.polynomialRing().getMonomial(m.exponents()),
					base.liftToInteger(s.getNumerator().coefficient(m)));
		}
		Map<Monomial, T> denominatorMap = new TreeMap<>();
		for (Monomial m : s.getDenominator().monomials()) {
			denominatorMap.put(extension.polynomialRing().getMonomial(m.exponents()),
					base.liftToInteger(s.getDenominator().coefficient(m)));
		}
		return extension.divide(extension.getEmbedding(extension.polynomialRing().getPolynomial(numeratorMap)),
				extension.getEmbedding(extension.polynomialRing().getPolynomial(denominatorMap)));
	}

	@Override
	public TExt<T> round(TExt<T> t, int accuracy) {
		Map<Monomial, T> numeratorMap = new TreeMap<>();
		for (Monomial m : t.getNumerator().monomials()) {
			numeratorMap.put(extension.polynomialRing().getMonomial(m.exponents()),
					base.round(t.getNumerator().coefficient(m), accuracy));
		}
		Map<Monomial, T> denominatorMap = new TreeMap<>();
		for (Monomial m : t.getDenominator().monomials()) {
			denominatorMap.put(extension.polynomialRing().getMonomial(m.exponents()),
					base.round(t.getDenominator().coefficient(m), accuracy));
		}
		return extension.divide(extension.getEmbedding(extension.polynomialRing().getPolynomial(numeratorMap)),
				extension.getEmbedding(extension.polynomialRing().getPolynomial(denominatorMap)));
	}

	@Override
	public int getAccuracy() {
		return base.getAccuracy();
	}

	@Override
	public TranscendentalDVRExtension<T, S> withAccuracy(int accuracy) {
		DiscreteValuationField<T, S> newBase = base.withAccuracy(accuracy);
		return new TranscendentalDVRExtension<>(newBase, extension.polynomialRing().getVariableNames());
	}

	@Override
	public ExtensionOfDiscreteValuationField<TExt<T>, TExt<S>, ?, ?, ?, ?, ?, ?, ?> getUniqueExtension(
			UnivariatePolynomial<TExt<T>> minimalPolynomial) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public OtherVersion<TExt<T>, ?, TExt<S>, ?> exact() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public OtherVersion<TExt<T>, ?, TExt<S>, ?> complete(int accuracy) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public TExt<T> getRandomInteger() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public DiscreteValuationRing<TExt<T>, TExt<S>> ringOfIntegers() {
		return localRing;
	}
}
