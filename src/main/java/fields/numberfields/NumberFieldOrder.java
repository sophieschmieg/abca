package fields.numberfields;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import fields.exceptions.InfinityException;
import fields.helper.AbstractAlgebra;
import fields.helper.AbstractIdeal;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.interfaces.Algebra;
import fields.interfaces.Ideal;
import fields.interfaces.Ring;
import fields.numberfields.NumberField.NFE;
import fields.numberfields.NumberFieldIntegers.NumberFieldIdeal;
import fields.vectors.FreeSubModule;
import fields.vectors.Matrix;
import fields.vectors.Vector;

public class NumberFieldOrder extends AbstractAlgebra<IntE, NFE> implements Algebra<IntE, NFE> {
	private NumberField field;
	private NumberFieldIntegers maximalOrder;
	private List<NFE> moduleBasis;
	private Matrix<IntE> latticedReducedToAsSubModule;
	private Matrix<IntE> asSubModuleToLatticeReduced;
	private FreeSubModule<IntE, NFE> asSubModule;

	NumberFieldOrder(NumberField field, List<NFE> generators) {
		this.field = field;
		this.maximalOrder = field.maximalOrder();
		List<NFE> moduleGenerators = new ArrayList<>();
		for (NFE generator : generators) {
			if (!field.isInteger(generator)) {
				throw new ArithmeticException("Not an integer generator!");
			}
			for (int i = 0; i < field.degree(); i++) {
				moduleGenerators.add(field.power(generator, i));
			}
		}
		this.asSubModule = new FreeSubModule<>(field.maximalOrder(), moduleGenerators);
		if (asSubModule.rank() != field.degree()) {
			throw new ArithmeticException("Not full rank!");
		}
		this.moduleBasis = maximalOrder.sublatticeReduction(asSubModule.getBasis());
		List<Vector<IntE>> latticeReducedInSubModule = new ArrayList<>();
		for (NFE latticeReduced : moduleBasis) {
			latticeReducedInSubModule.add(asSubModule.asVector(latticeReduced));
		}
		this.latticedReducedToAsSubModule = Matrix.fromColumns(latticeReducedInSubModule);
		this.asSubModuleToLatticeReduced = asSubModule.matrixAlgebra().inverse(latticedReducedToAsSubModule);
	}

	public IntE conductor() {
		return asSubModule.conductor();
	}

	public IntE discriminant() {
		return maximalOrder.discriminant(moduleBasis);
	}

	public boolean isMaximal() {
		return conductor().equals(Integers.z().one());
	}

	public boolean contains(NFE t) {
		return asSubModule.contains(t);
	}

	public boolean contains(NumberFieldOrder other) {
		return asSubModule.contains(other.asSubModule);
	}

	@Override
	public Ring<IntE> getRing() {
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
	public boolean isFree() {
		return true;
	}

	@Override
	public Ideal<IntE> annihilator() {
		return Integers.z().getZeroIdeal();
	}

	@Override
	public boolean isLinearIndependent(List<NFE> s) {
		return asSubModule.isLinearIndependent(s);
	}

	@Override
	public boolean isGeneratingModule(List<NFE> s) {
		return asSubModule.isGeneratingModule(s);
	}

	@Override
	public List<List<IntE>> nonTrivialCombinations(List<NFE> s) {
		return asSubModule.nonTrivialCombinations(s);
	}

	@Override
	public List<NFE> getModuleGenerators() {
		return moduleBasis;
	}
	
	@Override
	public Vector<IntE> asVector(NFE s) {
		return asSubModule.matrixAlgebra().multiply(asSubModuleToLatticeReduced, asSubModule.asVector(s));
	}

	@Override
	public Exactness exactness() {
		return field.exactness();
	}

	@Override
	public NFE getRandomElement() {
		return asSubModule.getRandomElement();
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
		return asSubModule.iterator();
	}

	@Override
	public NFE one() {
		return field.one();
	}

	@Override
	public BigInteger characteristic() {
		return BigInteger.ZERO;
	}

	@Override
	public NFE multiply(NFE t1, NFE t2) {
		return field.multiply(t1, t2);
	}

	@Override
	public boolean isUnit(NFE t) {
		return contains(field.inverse(t));
	}

	@Override
	public NFE inverse(NFE t) {
		NFE inverse = field.inverse(t);
		if (!contains(inverse)) {
			throw new ArithmeticException("Not a unit!");
		}
		return inverse;
	}

	@Override
	public boolean isCommutative() {
		return true;
	}

	@Override
	public boolean isIntegral() {
		return true;
	}

	@Override
	public boolean isReduced() {
		return true;
	}

	@Override
	public boolean isIrreducible() {
		return true;
	}

	@Override
	public boolean isZeroDivisor(NFE t) {
		return t.equals(zero());
	}

	@Override
	public boolean isEuclidean() {
		return false;
	}

	@Override
	public boolean isUniqueFactorizationDomain() {
		return false;
	}

	@Override
	public FactorizationResult<NFE, NFE> uniqueFactorization(NFE t) {
		throw new ArithmeticException("Not a UFD!");
	}

	@Override
	public boolean isPrincipalIdealDomain() {
		return false;
	}

	@Override
	public boolean isDedekindDomain() {
		return false;
	}

	@Override
	public boolean isDivisible(NFE dividend, NFE divisor) {
		if (divisor.equals(zero())) {
			return false;
		}
		return contains(field.divide(dividend, divisor));
	}

	@Override
	public QuotientAndRemainderResult<NFE> quotientAndRemainder(NFE dividend, NFE divisor) {
		if (divisor.equals(zero())) {
			return new QuotientAndRemainderResult<>(zero(), dividend);
		}
		NFE divided = field.divide(dividend, divisor);
		if (contains(divided)) {
			return new QuotientAndRemainderResult<>(divided, zero());
		}
		return new QuotientAndRemainderResult<>(zero(), dividend);
	}

	@Override
	public BigInteger euclidMeasure(NFE t) {
		return null;
	}

	@Override
	public NFE projectToUnit(NFE t) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Iterable<NFE> getUnits() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public int krullDimension() {
		return 1;
	}

	@Override
	public IdealResult<NFE, NumberFieldOrderIdeal> getIdealWithTransforms(List<NFE> generators) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Ideal<NFE> intersect(Ideal<NFE> t1, Ideal<NFE> t2) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Ideal<NFE> radical(Ideal<NFE> t) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public NFE getEmbedding(IntE t) {
		return maximalOrder.getEmbedding(t);
	}

	@Override
	public boolean isGeneratingAlgebra(List<NFE> s) {
		List<NFE> asModuleGenerators = new ArrayList<>();
		for (NFE generator : s) {
			for (int i = 0; i < field.degree(); i++) {
				asModuleGenerators.add(field.power(generator, i));
			}
		}
		return isGeneratingModule(asModuleGenerators);
	}

	@Override
	public List<NFE> getAlgebraGenerators() {
		return moduleBasis;
	}

	public class NumberFieldOrderIdeal extends AbstractIdeal<NFE> implements Ideal<NFE> {
		private FreeSubModule<IntE, NFE> asSubModule;
		private NumberFieldIdeal asNumberFieldIdeal;

		private NumberFieldOrderIdeal(List<NFE> generators) {
			super(NumberFieldOrder.this);
		}

		@Override
		public List<List<NFE>> nonTrivialCombinations(List<NFE> s) {
			// TODO Auto-generated method stub
			return null;
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
		public boolean isPrimary() {
			// TODO Auto-generated method stub
			return false;
		}

		@Override
		public boolean isPrime() {
			// TODO Auto-generated method stub
			return false;
		}

		@Override
		public boolean isMaximal() {
			// TODO Auto-generated method stub
			return false;
		}

		@Override
		public List<NFE> generators() {
			// TODO Auto-generated method stub
			return null;
		}

		@Override
		public List<NFE> generate(NFE t) {
			// TODO Auto-generated method stub
			return null;
		}

		@Override
		public NFE residue(NFE t) {
			// TODO Auto-generated method stub
			return null;
		}

		@Override
		public boolean contains(NFE t) {
			// TODO Auto-generated method stub
			return false;
		}
	}
}
