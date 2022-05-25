package fields.numberfields;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import fields.exceptions.InfinityException;
import fields.floatingpoint.Reals.Real;
import fields.helper.AbstractAlgebra;
import fields.helper.AbstractIdeal;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Rationals.Fraction;
import fields.interfaces.Algebra;
import fields.interfaces.Ideal;
import fields.interfaces.Lattice;
import fields.interfaces.RealInnerProductSpace;
import fields.interfaces.Ring;
import fields.numberfields.NumberField.NFE;
import fields.numberfields.NumberFieldIntegers.NumberFieldIdeal;
import fields.vectors.FreeModule;
import fields.vectors.FreeSubModule;
import fields.vectors.GenericPIDModule;
import fields.vectors.GenericPIDModule.Mod;
import fields.vectors.Matrix;
import fields.vectors.MatrixModule;
import fields.vectors.Vector;

public class NumberFieldOrder extends AbstractAlgebra<IntE, NFE>
		implements Algebra<IntE, NFE>, Lattice<NFE, Real, Vector<Real>> {
	private NumberField field;
	private NumberFieldIntegers maximalOrder;
	private List<NFE> moduleBasis;
	private Matrix<Real> generatorsAsMatrix;
	private Matrix<IntE> latticedReducedToAsSubModule;
	private Matrix<IntE> asSubModuleToLatticeReduced;
	private Matrix<Fraction> fromPowerBasisToModuleBasis;
	private FreeSubModule<IntE, NFE> asSubModule;
	private GenericPIDModule<IntE, NFE> mod;
	private GenericPIDModule<IntE, Vector<IntE>> minimalUnitPowerInOrder;
	private List<NFE> unitGroupGenerators;
	private OrderIdealGroup idealGroup;
	private PicardGroup picardGroup;

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
		this.mod = new GenericPIDModule<>(maximalOrder, asSubModule);
		this.moduleBasis = maximalOrder.sublatticeReduction(asSubModule.getBasis(), 1.0);
		this.asSubModule = new FreeSubModule<>(field.maximalOrder(), this.moduleBasis);
		List<Vector<IntE>> latticeReducedInSubModule = new ArrayList<>();
		List<Vector<Real>> generatorsAsVectors = new ArrayList<>();
		List<Vector<Fraction>> inPowerBasis = new ArrayList<>();
		for (NFE latticeReduced : moduleBasis) {
			latticeReducedInSubModule.add(asSubModule.asVector(latticeReduced));
			generatorsAsVectors.add(embedding(latticeReduced));
			inPowerBasis.add(field.asVector(latticeReduced));
		}
		Matrix<Fraction> moduleBasisToPowerBasis = Matrix.fromColumns(inPowerBasis);
		this.fromPowerBasisToModuleBasis = field.matrixAlgebra().inverse(moduleBasisToPowerBasis);
		this.latticedReducedToAsSubModule = Matrix.fromColumns(latticeReducedInSubModule);
		this.asSubModuleToLatticeReduced = asSubModule.matrixAlgebra().inverse(latticedReducedToAsSubModule);
		this.generatorsAsMatrix = Matrix.fromColumns(generatorsAsVectors);
	}

	@Override
	public String toString() {
		StringBuilder build = new StringBuilder();
		build.append("<");
		boolean first = true;
		for (NFE basis : moduleBasis) {
			if (first) {
				first = false;
			} else {
				build.append(", ");
			}
			build.append(basis);
		}
		build.append(">");
		return build.toString();
	}

	public NumberField numberField() {
		return field;
	}

	public NumberFieldIdeal conductorInMaximalOrder() {
		return maximalOrder.extend(mod.annihilator());
	}

	public NumberFieldOrderIdeal conductor() {
		return getIdeal(conductorInMaximalOrder().asSubModule().getBasis());
	}

	public IntE integerConductor() {
		return mod.annihilator().generators().get(0);
	}

	public IntE discriminant() {
		return maximalOrder.discriminant(moduleBasis);
	}

	public boolean isMaximal() {
		return integerConductor().equals(Integers.z().one());
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
	public List<Vector<IntE>> getSyzygies() {
		return Collections.emptyList();
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
	public List<Vector<IntE>> nonTrivialCombinations(List<NFE> s) {
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

	Vector<Fraction> asRationalVector(NFE s) {
		return field.matrixAlgebra().multiply(fromPowerBasisToModuleBasis, field.asVector(s));
	}

	public boolean isElement(NFE t) {
		return maximalOrder.isElement(t) && asSubModule.contains(t);
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
	public RealInnerProductSpace<Real, Vector<Real>> getVectorSpace() {
		return maximalOrder.getVectorSpace();
	}

	@Override
	public Vector<Real> embedding(NFE t) {
		return maximalOrder.embedding(t);
	}

	@Override
	public int rank() {
		return maximalOrder.rank();
	}

	@Override
	public Matrix<Real> generatorsAsMatrix() {
		return generatorsAsMatrix;
	}

	public OrderIdealGroup idealGroup() {
		if (idealGroup == null) {
			idealGroup = new OrderIdealGroup(this);
		}
		return idealGroup;
	}

	public PicardGroup picardGroup() {
		if (picardGroup == null) {
			picardGroup = new PicardGroup(this);
		}
		return picardGroup;
	}

	public BigInteger classNumber() {
		return picardGroup().getNumberOfElements();
		// BigInteger cl = maximalOrder.numberField().classNumber();
//		BigInteger unitsMaximalOrderModConductor = conductorInMaximalOrder().modOut().getNumberOfUnits();
//		// BigInteger unitsOrderModConductor = conductor().modOut().getNumberOfUnits();
//		GenericPIDModule<IntE, Vector<IntE>> maximalUnits = maximalOrder.unitGroupAsIntModule();
//		return null;
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

	public List<NFE> unitGroupGenerators() {
		if (unitGroupGenerators == null) {
			Integers z = Integers.z();
			List<NFE> maximalOrderUnitGroupGenerators = maximalOrder.unitGroupGenerators();
			FreeModule<IntE> free = new FreeModule<>(z, maximalOrderUnitGroupGenerators.size());
			List<Vector<IntE>> minimumPower = new ArrayList<>();
			unitGroupGenerators = new ArrayList<>();
			for (int i = 0; i < maximalOrderUnitGroupGenerators.size(); i++) {
				NFE unit = maximalOrderUnitGroupGenerators.get(i);
				int power = 1;
				NFE elementUnit = unit;
				while (!isElement(elementUnit)) {
					power++;
					elementUnit = maximalOrder.multiply(unit, elementUnit);
				}
				unitGroupGenerators.add(elementUnit);
				minimumPower.add(free.scalarMultiply(power, free.getUnitVector(i + 1)));
			}
			minimalUnitPowerInOrder = new GenericPIDModule<>(free, new FreeSubModule<>(free, minimumPower));
		}
		return unitGroupGenerators;
	}

	Iterator<NFE> unitClasses() {
		unitGroupGenerators();
		return new Iterator<>() {
			private Iterator<Mod<Vector<IntE>>> it = minimalUnitPowerInOrder.iterator();

			@Override
			public boolean hasNext() {
				return it.hasNext();
			}

			@Override
			public NFE next() {
				return maximalOrder.fromUnitGeneratorVector(minimalUnitPowerInOrder.lift(it.next()));
			}
		};
	}

	@Override
	public int krullDimension() {
		return 1;
	}

	private Matrix<IntE> asMatrix(NFE t) {
		List<Vector<IntE>> columns = new ArrayList<>();
		for (NFE integral : getModuleGenerators()) {
			columns.add(asVector(multiply(t, integral)));
		}
		return Matrix.fromColumns(columns);
	}

	private Matrix<IntE> asIntegerMatrix(Matrix<NFE> m) {
		List<List<Matrix<IntE>>> result = new ArrayList<>();
		for (int i = 0; i < m.rows(); i++) {
			List<Matrix<IntE>> row = new ArrayList<>();
			for (int j = 0; j < m.columns(); j++) {
				row.add(asMatrix(m.entry(i + 1, j + 1)));
			}
		}
		return Matrix.fromBlockMatrix(result);
	}

	private Vector<IntE> asIntegerVector(Vector<NFE> b) {
		List<IntE> result = new ArrayList<>();
		for (NFE t : b.asList()) {
			result.addAll(asVector(t).asList());
		}
		return new Vector<IntE>(result);
	}

	private Vector<NFE> fromIntegerVector(Vector<IntE> b) {
		List<NFE> result = new ArrayList<>();
		for (int i = 0; i < b.dimension() / field.degree(); i++) {
			result.add(fromVector(new Vector<>(b.asList().subList(i * field.degree(), (i + 1) * field.degree()))));
		}
		return new Vector<NFE>(result);
	}

	@Override
	public boolean isSubModuleMember(MatrixModule<NFE> module, Matrix<NFE> m, Vector<NFE> b) {
		Integers z = Integers.z();
		return z.isSubModuleMember(new MatrixModule<>(z, field.degree() * m.rows(), field.degree() * m.columns()),
				asIntegerMatrix(m), asIntegerVector(b));
	}

	@Override
	public Vector<NFE> asSubModuleMember(MatrixModule<NFE> module, Matrix<NFE> m, Vector<NFE> b) {
		Integers z = Integers.z();
		return fromIntegerVector(
				z.asSubModuleMember(new MatrixModule<>(z, field.degree() * m.rows(), field.degree() * m.columns()),
						asIntegerMatrix(m), asIntegerVector(b)));
	}

	@Override
	public List<Vector<NFE>> syzygyProblem(MatrixModule<NFE> module, Matrix<NFE> m) {
		Integers z = Integers.z();
		List<Vector<IntE>> integerSyzygies = z.syzygyProblem(
				new MatrixModule<>(z, field.degree() * m.rows(), field.degree() * m.columns()), asIntegerMatrix(m));
		List<Vector<NFE>> syzygies = new ArrayList<>();
		for (Vector<IntE> integerSyzygy : integerSyzygies) {
			syzygies.add(fromIntegerVector(integerSyzygy));
		}
		return syzygies;
	}

	@Override
	public List<Vector<NFE>> simplifySubModuleGenerators(MatrixModule<NFE> module, Matrix<NFE> m) {
		List<Vector<NFE>> result = new ArrayList<>();
		for (int i = 0; i < m.columns(); i++) {
			Vector<NFE> column = m.column(i + 1);
			if (!column.equals(module.codomain().zero())
					&& (result.isEmpty() || !isSubModuleMember(Matrix.fromColumns(result), column))) {
				result.add(column);
			}
		}
		return result;
	}

	@Override
	public IdealResult<NFE, NumberFieldOrderIdeal> getIdealWithTransforms(List<NFE> generators) {
		NumberFieldOrderIdeal ideal = getIdeal(generators);
		List<NFE> integralBasis = getModuleGenerators();
		List<Vector<IntE>> columns = new ArrayList<>();
		for (NFE generator : generators) {
			for (NFE basisElement : integralBasis) {
				NFE multiplied = multiply(basisElement, generator);
				columns.add(asVector(multiplied));
			}
		}
		Matrix<IntE> matrix = Matrix.fromColumns(columns);
		MatrixModule<IntE> mm = new MatrixModule<>(Integers.z(), integralBasis.size(), columns.size());
		List<List<NFE>> transforms = new ArrayList<>();
		for (NFE basisElement : ideal.generators()) {
			Vector<IntE> rhs = asVector(basisElement);
			List<IntE> solution = mm.solve(matrix, rhs).asList();
			List<NFE> transform = new ArrayList<>();
			for (int i = 0; i < generators.size(); i++) {
				transform.add(fromVector(
						new Vector<>(solution.subList(i * integralBasis.size(), (i + 1) * integralBasis.size()))));
			}
			transforms.add(transform);
		}
		Matrix<NFE> row = Matrix.fromRows(Collections.singletonList(new Vector<>(generators)));
		return new IdealResult<>(transforms, generators, ideal, syzygyProblem(row.getModule(this), row));
	}

	@Override
	public NumberFieldOrderIdeal getIdeal(List<NFE> generators) {
		return new NumberFieldOrderIdeal(generators);
	}

	@Override
	public NumberFieldOrderIdeal getIdeal(NFE... generators) {
		return getIdeal(Arrays.asList(generators));
	}

	@Override
	public NumberFieldOrderIdeal getUnitIdeal() {
		return getIdeal(one());
	}

	@Override
	public NumberFieldOrderIdeal getZeroIdeal() {
		return getIdeal(zero());
	}

	@Override
	public NumberFieldOrderIdeal intersect(Ideal<NFE> t1, Ideal<NFE> t2) {
		NumberFieldOrderIdeal ideal1 = (NumberFieldOrderIdeal) t1;
		NumberFieldOrderIdeal ideal2 = (NumberFieldOrderIdeal) t2;
		return getIdeal(ideal1.asSubModule.intersection(ideal2.asSubModule).getBasis());
	}

	@Override
	public NumberFieldOrderIdeal radical(Ideal<NFE> t) {
		// TODO Auto-generated method stub
		return null;
	}

	public boolean extensibleToMaximalOrder(NumberFieldOrderIdeal t) {
		return coprime(t, conductor());
	}

	public NumberFieldIdeal extendToMaximalOrder(NumberFieldOrderIdeal t) {
		return maximalOrder.getIdeal(t.generators());
	}

	public boolean restrictableFromMaximalOrder(NumberFieldIdeal t) {
		return maximalOrder.coprime(t, conductorInMaximalOrder());
	}

	public NumberFieldOrderIdeal restrictFromMaximalOrder(NumberFieldIdeal t) {
		return getIdeal(asSubModule.intersection(t.asSubModule()).getBasis());
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
		private List<NFE> generators;
		private FreeSubModule<IntE, NFE> asSubModule;
		private Matrix<IntE> generatorsToAsSubModule;
		private Matrix<IntE> asSubModuleToGenerators;
		private GenericPIDModule<IntE, NFE> mod;
		private boolean principalComputed;
		private NFE principalGenerator;

		private NumberFieldOrderIdeal(List<NFE> generators) {
			super(NumberFieldOrder.this);
			NumberFieldOrder order = NumberFieldOrder.this;
			List<NFE> multiplied = new ArrayList<>();
			for (NFE b : order.getModuleGenerators()) {
				for (NFE generator : generators) {
					multiplied.add(multiply(b, generator));
				}
			}
			this.asSubModule = new FreeSubModule<>(order, multiplied);
			if (asSubModule.rank() == 0) {
				this.generators = Collections.emptyList();
				return;
			}
			this.generators = order.getVectorSpace().latticeReduction(asSubModule.getBasis(), order, 1.0);
			this.asSubModule = new FreeSubModule<>(order, this.generators);
			List<Vector<IntE>> generatorList = new ArrayList<>();
			for (NFE generator : this.generators) {
				generatorList.add(asSubModule.asVector(generator));
			}
			this.generatorsToAsSubModule = Matrix.fromColumns(generatorList);
			this.asSubModuleToGenerators = asSubModule.matrixAlgebra().inverse(generatorsToAsSubModule);
		}

		private GenericPIDModule<IntE, NFE> mod() {
			if (mod == null) {
				mod = new GenericPIDModule<>(NumberFieldOrder.this, asSubModule);
			}
			return mod;
		}

		public FreeSubModule<IntE, NFE> asSubModule() {
			return asSubModule;
		}

//		@Override
//		public List<List<NFE>> nonTrivialCombinations(List<NFE> s) {
//			List<Vector<IntE>> asVectorList = new ArrayList<>();
//			for (NFE e : s) {
//				for (NFE basisVector : NumberFieldOrder.this.getModuleGenerators()) {
//					asVectorList.add(NumberFieldOrder.this.asVector(multiply(e, basisVector)));
//				}
//			}
//			List<List<IntE>> nonTrivialCombinations = new FreeModule<>(Integers.z(), maximalOrder.rank())
//					.nonTrivialCombinations(asVectorList);
//			List<List<NFE>> result = new ArrayList<>();
//			int degree = field.degree();
//			for (List<IntE> combination : nonTrivialCombinations) {
//				List<NFE> row = new ArrayList<>();
//				for (int i = 0; i < s.size(); i++) {
//					row.add(NumberFieldOrder.this
//							.fromVector(new Vector<>(combination.subList(i * degree, (i + 1) * degree))));
//				}
//				boolean found = false;
//				for (List<NFE> prevRow : result) {
//					NFE ratio = null;
//					boolean notFound = false;
//					for (int i = 0; i < s.size(); i++) {
//						QuotientAndRemainderResult<NFE> qr = quotientAndRemainder(row.get(i), prevRow.get(i));
//						if (!qr.getRemainder().equals(zero())) {
//							notFound = true;
//							break;
//						}
//						if (ratio == null) {
//							ratio = qr.getQuotient();
//						}
//						if (!ratio.equals(qr.getQuotient())) {
//							notFound = true;
//							break;
//						}
//					}
//					found = !notFound;
//					if (found) {
//						break;
//					}
//				}
//				if (!found) {
//					result.add(row);
//				}
//			}
//			return result;
//		}

		@Override
		public boolean isFinite() {
			return false;
		}

		@Override
		public BigInteger getNumberOfElements() throws InfinityException {
			throw new InfinityException();
		}

		@Override
		public boolean isPrincipal() {
			if (!principalComputed) {
				principalComputed = true;
				if (generators.isEmpty()) {
					principalGenerator = zero();
					return true;
				}
				NumberFieldIdeal extended = extendToMaximalOrder(this);
				if (!extended.isPrincipal()) {
					return false;
				}
				NFE generator = extended.principalGenerator();
				Iterator<NFE> it = unitClasses();
				while (it.hasNext()) {
					NFE unit = it.next();
					NFE adjusted = maximalOrder.multiply(unit, generator);
					if (!isElement(adjusted) || !contains(adjusted)) {
						continue;
					}
					boolean principal = true;
					for (NFE g : generators) {
						if (!isDivisible(g, adjusted)) {
							principal = false;
							break;
						}
					}
					if (principal) {
						principalGenerator = adjusted;
						return true;
					}
				}
			}
			return principalGenerator != null;
		}

		public NFE principalGenerator() {
			isPrincipal();
			return principalGenerator;
		}

		@Override
		public boolean isPrimary() {
			if (generators.isEmpty()) {
				return true;
			}
			// TODO Auto-generated method stub
			return false;
		}

		@Override
		public boolean isPrime() {
			if (generators.isEmpty()) {
				return true;
			}
			if (extensibleToMaximalOrder(this)) {
				return extendToMaximalOrder(this).isPrime();
			}
			// TODO Auto-generated method stub
			return false;
		}

		@Override
		public boolean isMaximal() {
			if (generators.isEmpty()) {
				return false;
			}
			// TODO Auto-generated method stub
			return false;
		}

		@Override
		public List<NFE> generators() {
			return generators;
		}

		@Override
		public List<NFE> generate(NFE t) {
			if (generators.isEmpty()) {
				return Collections.emptyList();
			}
			Vector<IntE> result = asSubModule.matrixAlgebra().multiply(asSubModuleToGenerators,
					asSubModule.asVector(t));
			List<NFE> asNFE = new ArrayList<>();
			for (IntE r : result.asList()) {
				asNFE.add(NumberFieldOrder.this.getInteger(r));
			}
			return asNFE;
		}

		@Override
		public NFE residue(NFE t) {
			if (generators.isEmpty()) {
				return t;
			}
			return mod().lift(mod().reduce(t));
		}

		@Override
		public boolean contains(NFE t) {
			if (generators.isEmpty()) {
				return t.equals(zero());
			}
			return asSubModule.contains(t);
		}
	}
}
