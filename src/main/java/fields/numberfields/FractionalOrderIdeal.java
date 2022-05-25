package fields.numberfields;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import fields.exceptions.InfinityException;
import fields.floatingpoint.FiniteRealVectorSpace;
import fields.floatingpoint.Reals.Real;
import fields.helper.AbstractModule;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Rationals.Fraction;
import fields.interfaces.Element;
import fields.interfaces.Ideal;
import fields.interfaces.Lattice;
import fields.numberfields.NumberField.NFE;
import fields.numberfields.NumberFieldOrder.NumberFieldOrderIdeal;
import fields.vectors.FreeModule;
import fields.vectors.FreeSubModule;
import fields.vectors.Matrix;
import fields.vectors.Vector;

public class FractionalOrderIdeal extends AbstractModule<IntE, NFE>
		implements Element<FractionalOrderIdeal>, Lattice<NFE, Real, Vector<Real>> {
	private NumberField field;
	private NumberFieldOrder order;
	private NumberFieldOrderIdeal numerator;
	private IntE denominator;
	private List<NFE> basis;
	private Matrix<Fraction> toBasisBaseChange;
	private Matrix<Fraction> fromBasisBaseChange;
	private Matrix<Real> generatorsAsMatrix;

	public FractionalOrderIdeal(NumberFieldOrder order, NumberFieldOrderIdeal numerator) {
		this(order, numerator.generators());
	}

	public FractionalOrderIdeal(NumberFieldOrder order, NumberFieldOrderIdeal numerator, IntE denominator) {
		this(order, numerator, order.getIdeal(Collections.singletonList(order.getInteger(denominator))));
	}

	private static List<NFE> generators(NumberFieldOrder order, NumberFieldOrderIdeal numerator,
			NumberFieldOrderIdeal denominator) {
		Integers z = Integers.z();
		List<NFE> inverted = new ArrayList<>();
		IntE integerDenominator = z.one();
		for (NFE denominatorGenerator : denominator.generators()) {
			NFE invertedDenominator = order.numberField().inverse(denominatorGenerator);
			inverted.add(invertedDenominator);
			Vector<Fraction> asFractions = order.asRationalVector(invertedDenominator);
			for (Fraction c : asFractions.asList()) {
				integerDenominator = z.lcm(c.getDenominator(), integerDenominator);
			}
		}
		NumberFieldOrderIdeal clearedInverseDenominator = order.getUnitIdeal();
		for (NFE invertedGenerator : inverted) {
			clearedInverseDenominator = order.intersect(
					order.getIdeal(order.numberField().multiply(integerDenominator, invertedGenerator)),
					clearedInverseDenominator);
		}
		List<NFE> generators = new ArrayList<>();
		for (NFE generator : order.multiply(numerator, clearedInverseDenominator).generators()) {
			generators.add(order.numberField().divide(generator, order.numberField().getInteger(integerDenominator)));
		}
		return generators;
	}

	public FractionalOrderIdeal(NumberFieldOrder order, NumberFieldOrderIdeal numerator,
			NumberFieldOrderIdeal denominator) {
		this(order, generators(order, numerator, denominator));
	}

	public FractionalOrderIdeal(NumberFieldOrder order, List<NFE> generators) {
		Integers z = Integers.z();
		this.order = order;
		this.field = order.numberField();
		denominator = z.one();
		List<NFE> multipliedGenerators = new ArrayList<>();
		for (NFE orderBasis : order.getModuleGenerators()) {
			for (NFE generator : generators) {
				NFE multiplied = order.multiply(orderBasis, generator);
				multipliedGenerators.add(multiplied);
				Vector<Fraction> asFractions = order.asRationalVector(multiplied);
				for (Fraction c : asFractions.asList()) {
					denominator = z.lcm(c.getDenominator(), denominator);
				}
			}
		}
		List<NFE> clearedGenerators = new ArrayList<>();
		for (NFE generator : multipliedGenerators) {
			clearedGenerators.add(field.multiply(denominator, generator));
		}
		clearedGenerators = new FreeSubModule<>(order, clearedGenerators).getBasis();
		if (clearedGenerators.isEmpty()) {
			this.numerator = order.getZeroIdeal();
			this.denominator = z.one();
			this.basis = Collections.emptyList();
			return;
		}
		clearedGenerators = field.minkowskiEmbeddingSpace().latticeReduction(clearedGenerators, this, 1.0);
		this.numerator = order.getIdeal(clearedGenerators);
		this.basis = new ArrayList<>();
		List<Vector<Fraction>> basisAsVectors = new ArrayList<>();
		for (NFE basisVector : clearedGenerators) {
			NFE divided = field.divide(basisVector, field.getInteger(denominator));
			basis.add(divided);
			basisAsVectors.add(field.asVector(divided));
		}
		this.fromBasisBaseChange = Matrix.fromColumns(basisAsVectors);
		this.toBasisBaseChange = field.matrixAlgebra().inverse(fromBasisBaseChange);
	}

	public static FractionalOrderIdeal principalIdeal(NumberFieldOrder order, NFE generator) {
		return new FractionalOrderIdeal(order, Collections.singletonList(generator));
	}

	@Override
	public String toString() {
		return getNumerator() + "/" + getDenominator();
	}

	public NumberFieldOrder getOrder() {
		return order;
	}

	public NumberFieldOrderIdeal getNumerator() {
		return numerator;
	}

	public IntE getDenominator() {
		return denominator;
	}

	public boolean equals(Object o) {
		if (!(o instanceof FractionalOrderIdeal)) {
			return false;
		}
		return compareTo((FractionalOrderIdeal) o) == 0;
	}

	@Override
	public int compareTo(FractionalOrderIdeal o) {
		int cmp = denominator.compareTo(o.denominator);
		if (cmp != 0) {
			return cmp;
		}
		return numerator.compareTo(o.numerator);
	}

	public boolean isPrincipal() {
		return numerator.isPrincipal();
	}

	public boolean isInteger() {
		return denominator.equals(Integers.z().one());
	}

	@Override
	public Integers getRing() {
		return Integers.z();
	}

	@Override
	public NFE zero() {
		return field.zero();
	}

	@Override
	public NFE add(NFE s1, NFE s2) {
		return field.add(s1, s2);
	}

	@Override
	public NFE negative(NFE s) {
		return field.negative(s);
	}

	@Override
	public NFE scalarMultiply(IntE t, NFE s) {
		return field.multiply(t, s);
	}

	@Override
	public boolean isFree() {
		return true;
	}

	@Override
	public Ideal<IntE> annihilator() {
		return Integers.z().getZeroIdeal();
	}
	
	@Override
	public List<Vector<IntE>> getSyzygies() {
		return Collections.emptyList();
	}

	private List<Vector<IntE>> asVectors(List<NFE> s) {
		List<Vector<IntE>> result = new ArrayList<>();
		for (NFE t : s) {
			result.add(asVector(t));
		}
		return result;
	}

	@Override
	public boolean isLinearIndependent(List<NFE> s) {
		return new FreeModule<>(Integers.z(), field.degree()).isLinearIndependent(asVectors(s));
	}

	@Override
	public boolean isGeneratingModule(List<NFE> s) {
		return new FreeModule<>(Integers.z(), field.degree()).isGeneratingModule(asVectors(s));
	}

	@Override
	public List<Vector<IntE>> nonTrivialCombinations(List<NFE> s) {
		return new FreeModule<>(Integers.z(), field.degree()).nonTrivialCombinations(asVectors(s));
	}

	@Override
	public List<NFE> getModuleGenerators() {
		return basis;
	}

	@Override
	public Matrix<Real> generatorsAsMatrix() {
		if (generatorsAsMatrix == null) {
			List<Vector<Real>> asVectors = new ArrayList<>();
			for (NFE generator : getModuleGenerators()) {
				asVectors.add(embedding(generator));
			}
			generatorsAsMatrix = Matrix.fromColumns(asVectors);
		}
		return generatorsAsMatrix;
	}

	@Override
	public Vector<IntE> asVector(NFE s) {
		Vector<Fraction> inBasis = field.matrixAlgebra().multiply(toBasisBaseChange, field.asVector(s));
		List<IntE> result = new ArrayList<>();
		for (Fraction c : inBasis.asList()) {
			result.add(c.asInteger());
		}
		return new Vector<>(result);
	}

	@Override
	public Exactness exactness() {
		return Exactness.EXACT;
	}

	@Override
	public NFE getRandomElement() {
		return fromVector(new FreeModule<IntE>(Integers.z(), field.degree()).getRandomElement());
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
	public Iterator<NFE> iterator() {
		throw new InfinityException();
	}

	@Override
	public FiniteRealVectorSpace getVectorSpace() {
		return field.minkowskiEmbeddingSpace();
	}

	@Override
	public Vector<Real> embedding(NFE t) {
		return field.minkowskiEmbedding(t);
	}

	@Override
	public int rank() {
		return field.degree();
	}
}
