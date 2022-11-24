package fields.numberfields;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.TreeSet;

import fields.exceptions.InfinityException;
import fields.finitefields.FiniteField;
import fields.finitefields.FiniteField.FFE;
import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.floatingpoint.Complex;
import fields.floatingpoint.Complex.ComplexNumber;
import fields.floatingpoint.FiniteRealVectorSpace;
import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.helper.AbstractAlgebra;
import fields.helper.AbstractIdeal;
import fields.helper.FieldEmbedding;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Integers.IntegerIdeal;
import fields.integers.LocalizedFractions;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.Algebra;
import fields.interfaces.DedekindRing;
import fields.interfaces.DedekindRingExtension;
import fields.interfaces.DiscreteValuationRing;
import fields.interfaces.DiscreteValuationRing.OkutsuType;
import fields.interfaces.DiscreteValuationRing.TheMontesResult;
import fields.interfaces.Ideal;
import fields.interfaces.Lattice;
import fields.interfaces.MathMap;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.local.LocalRingExtension;
import fields.local.Value;
import fields.numberfields.ModuloNumberFieldIdeal.ModNFE;
import fields.numberfields.NumberField.NFE;
import fields.vectors.DualVectorSpace;
import fields.vectors.DualVectorSpace.Dual;
import fields.vectors.FreeModule;
import fields.vectors.FreeSubModule;
import fields.vectors.GenericPIDModule;
import fields.vectors.Matrix;
import fields.vectors.MatrixModule;
import fields.vectors.Polytope;
import fields.vectors.RealLattice;
import fields.vectors.RealNumberFieldIntegerLattice;
import fields.vectors.Vector;
import util.FunctionMathMap;
import util.Identity;
import util.MiscAlgorithms;
import util.SingletonSortedMap;

public class NumberFieldIntegers extends AbstractAlgebra<IntE, NFE> implements Algebra<IntE, NFE>,
		DedekindRingExtension<Fraction, IntE, PFE, NFE, NFE, PFE, FFE, FiniteField, LocalizedNumberField, NumberField>,
		Lattice<NFE, Real, Vector<Real>> {
	private NumberField field;
	private UnivariatePolynomial<IntE> minimalPolynomial;
	private IntE minimalPolynomialDenominator;
	private List<NFE> integralBasis;
	private Matrix<Fraction> toIntegralBasis;
	private Matrix<Real> generatorsAsMatrix;
	private Matrix<Fraction> fromIntegralBasis;
	private Set<NFE> units;
	private SortedMap<IntE, NFE> primitiveRootsOfUnity;
	private RealLattice unitLattice;
	private List<NFE> unitIdealMultipliers;
	private Map<IntE, List<NumberFieldIdeal>> idealsOverPrime;
	private Map<FactorizationResult<Ideal<NFE>, Ideal<NFE>>, NumberFieldIdeal> idealsByFactorization;
	private NumberFieldIdeal unitIdeal;
	private NumberFieldIdeal zeroIdeal;
	private FreeModule<IntE> asFreeModule;
	private Map<Ideal<NFE>, LocalizedNumberField> localizations;

	NumberFieldIntegers(NumberField field) {
		this.field = field;
		this.idealsOverPrime = new TreeMap<>();
		this.idealsByFactorization = new TreeMap<>();
		this.asFreeModule = new FreeModule<>(Integers.z(), field.degree());
		this.localizations = new TreeMap<>();
	}

	@Override
	public String toString() {
		return "O_" + field.toString();
	}

	@Override
	public Exactness exactness() {
		return field.exactness();
	}

	@Override
	public NFE getRandomElement() {
		return fromVector(asFreeModule.getRandomElement());
	}

	@Override
	public boolean isCommutative() {
		return true;
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
	public NFE zero() {
		return field.zero();
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
	public NFE add(NFE t1, NFE t2) {
		return field.add(t1, t2);
	}

	@Override
	public NFE negative(NFE t) {
		return field.negative(t);
	}

	@Override
	public NFE multiply(NFE t1, NFE t2) {
		return field.multiply(t1, t2);
	}

	@Override
	public boolean isUnit(NFE t) {
		return Integers.z().isUnit(field.norm(t).getNumerator());
	}

	@Override
	public NFE inverse(NFE t) {
		if (!isUnit(t)) {
			throw new ArithmeticException(t + " is not invertible");
		}
		return field.inverse(t);
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
		return false;
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
		return true;
	}

	@Override
	public boolean isDivisible(NFE dividend, NFE divisor) {
		return isElement(field.divide(dividend, divisor));
	}

	@Override
	public QuotientAndRemainderResult<NFE> quotientAndRemainder(NFE dividend, NFE divisor) {
		Ideal<NFE> divisorIdeal = getIdeal(Collections.singletonList(divisor));
		NFE residue = divisorIdeal.residue(dividend);
		return new QuotientAndRemainderResult<>(field.divide(subtract(dividend, residue), divisor), residue);
	}

	@Override
	public BigInteger euclidMeasure(NFE t) {
		throw new ArithmeticException("Not a Euclidean ring!");
	}

	@Override
	public NFE projectToUnit(NFE t) {
		NumberFieldIdeal ideal = getIdeal(Collections.singletonList(t));
		NFE maxGenerator = ideal.principalGenerator();
		int max = primitiveRootsOfUnity().lastKey().intValueExact();
		NFE primitiveRootOfUnity = primitiveRootsOfUnity().get(Integers.z().getInteger(max));
		NFE generator = maxGenerator;
		for (int i = 0; i < max - 1; i++) {
			generator = multiply(primitiveRootOfUnity, generator);
			if (generator.compareTo(maxGenerator) > 0) {
				maxGenerator = generator;
			}
		}
		return divideChecked(t, maxGenerator);
	}

	@Override
	public Iterable<NFE> getUnits() {
		if (field.degree() > 2) {
			throw new UnsupportedOperationException("Not implemented!");
		}
		if (field.degree() == 1) {
			List<NFE> units = new ArrayList<>();
			units.add(getInteger(-1));
			units.add(getInteger(1));
			return units;
		}
		Integers z = Integers.z();
		if (field.discriminant().compareTo(z.zero()) < 0) {
			List<NFE> units = new ArrayList<>();
			units.add(getInteger(-1));
			units.add(getInteger(1));
			if (field.discriminant().equals(z.negative(z.getInteger(4)))) {
				units.addAll(field.sqrt(field.getInteger(-1)).keySet());
			} else if (field.discriminant().equals(z.negative(z.getInteger(3)))) {
				UnivariatePolynomialRing<NFE> polynomials = field.getUnivariatePolynomialRing();
				UnivariatePolynomial<NFE> eisenstein = polynomials.getPolynomial(field.one(), field.one(), field.one());
				Set<NFE> zeta3 = field.roots(eisenstein).keySet();
				units.addAll(zeta3);
				for (NFE unit : zeta3) {
					units.add(field.negative(unit));
				}
			}
			return units;
		}
		IntE discriminant = field.discriminant();
		boolean mod23;
		if (discriminant.getValue().mod(BigInteger.valueOf(4)).equals(BigInteger.ZERO)) {
			discriminant = Integers.z().divideChecked(discriminant, z.getInteger(4));
			mod23 = true;
		} else {
			mod23 = false;
		}
		List<Vector<IntE>> pells = MiscAlgorithms.pellsEquation(discriminant, z.getInteger(mod23 ? -1 : -4), false);
		pells.addAll(MiscAlgorithms.pellsEquation(discriminant, z.getInteger(mod23 ? 1 : 4), false));
		Vector<IntE> pell = pells.get(0);
		NFE alpha = field.sqrt(field.getInteger(discriminant)).keySet().iterator().next();
		final NFE fundamental;
		if (mod23) {
			fundamental = add(getInteger(pell.get(1)), multiply(pell.get(2), alpha));
		} else {
			fundamental = divide(add(getInteger(pell.get(1)), multiply(pell.get(2), alpha)), getInteger(2));
		}
		return new Iterable<NumberField.NFE>() {

			@Override
			public Iterator<NFE> iterator() {
				return new Iterator<>() {
					Iterator<IntE> it = z.iterator();
					IntE current = null;

					@Override
					public boolean hasNext() {
						return true;
					}

					@Override
					public NFE next() {
						if (current == null) {
							current = it.next();
							return power(fundamental, current.intValueExact());
						}
						NFE result = negative(power(fundamental, current.intValueExact()));
						current = null;
						return result;
					}
				};
			}
		};
	}

	public SortedMap<IntE, NFE> primitiveRootsOfUnity() {
		if (primitiveRootsOfUnity == null) {
			Integers z = Integers.z();
			UnivariatePolynomialRing<NFE> polynomials = field.getUnivariatePolynomialRing();
			primitiveRootsOfUnity = new TreeMap<>();
			primitiveRootsOfUnity.put(z.one(), one());
			List<List<IntE>> primePowers = new ArrayList<>();
			for (IntE prime : z.setOfPrimes()) {
				if (prime.compareTo(z.getInteger(field.degree() + 1)) > 0) {
					break;
				}
				UnivariatePolynomial<NFE> unity = polynomials.toUnivariate(polynomials.divideChecked(
						polynomials.subtract(polynomials.getVarPower(prime.intValueExact()), polynomials.one()),
						polynomials.subtract(polynomials.getVar(), polynomials.one())));
				Map<NFE, Integer> roots = field.roots(unity);
				if (roots.isEmpty()) {
					continue;
				}
				NFE rootOfUnity = roots.keySet().iterator().next();
				primitiveRootsOfUnity.put(prime, rootOfUnity);
				List<IntE> powers = new ArrayList<>();
				powers.add(z.one());
				powers.add(prime);
				for (int i = 2; z.multiply(prime.intValueExact() - 1, z.power(prime, i - 1)).intValueExact() <= field
						.degree(); i++) {
					roots = field.roots(rootOfUnity, prime.intValueExact());
					if (roots.isEmpty()) {
						break;
					}
					IntE primePower = z.power(prime, i);
					rootOfUnity = roots.keySet().iterator().next();
					primitiveRootsOfUnity.put(primePower, rootOfUnity);
					powers.add(primePower);
				}
				primePowers.add(powers);
			}
			List<List<IntE>> factors = MiscAlgorithms.crossProduct(primePowers);
			for (List<IntE> factorList : factors) {
				IntE factor = z.one();
				NFE rootOfUnity = one();
				for (IntE primePower : factorList) {
					factor = z.multiply(primePower, factor);
					rootOfUnity = multiply(primitiveRootsOfUnity.get(primePower), rootOfUnity);
				}
				primitiveRootsOfUnity.put(factor, rootOfUnity);
			}
		}
		return primitiveRootsOfUnity;
	}

	public Set<NFE> freeUnitGroupGenerators() {
		if (units == null) {
			int numberOfUnits = field.realEmbeddings().size() + field.complexEmbeddings().size() - 1;
			if (numberOfUnits == 0) {
				units = Collections.emptySet();
			} else {
				idealEquivalenceClass(getUnitIdeal());
			}
		}
		return units;
	}

	public List<NFE> unitGroupGenerators() {
		List<NFE> result = new ArrayList<>();
		SortedMap<IntE, NFE> rootsOfUnity = primitiveRootsOfUnity();
		result.add(rootsOfUnity.get(rootsOfUnity.lastKey()));
		result.addAll(freeUnitGroupGenerators());
		return result;
	}

	public Vector<IntE> asUnitGeneratorVector(NFE t) {
		if (!isUnit(t)) {
			throw new ArithmeticException("Not a unit!");
		}
		Integers z = Integers.z();
		int numberOfUnitGenerators = unitGroupGenerators().size();
		NFE rootOfUnity;
		Vector<IntE> free = new Vector<>();
		if (numberOfUnitGenerators > 1) {
			Vector<Real> logs = field.logRepresentation(t);
			free = unitLattice.asVector(logs);
			List<IntE> freeList = new ArrayList<>();
			freeList.add(z.one());
			freeList.addAll(free.asList());
			NFE fromFree = fromUnitGeneratorVector(new Vector<>(freeList));
			rootOfUnity = divide(t, fromFree);
		} else {
			rootOfUnity = t;
		}
		int max = primitiveRootsOfUnity().lastKey().intValueExact();
		NFE primitiveRootInverted = inverse(primitiveRootsOfUnity().get(primitiveRootsOfUnity().lastKey()));
		for (int i = 0; i < max; i++) {
			if (rootOfUnity.equals(one())) {
				List<IntE> result = new ArrayList<>();
				result.add(z.getInteger(i));
				result.addAll(free.asList());
				return new Vector<>(result);
			}
			rootOfUnity = multiply(rootOfUnity, primitiveRootInverted);
		}
		throw new ArithmeticException("Did not get a root of unity!");
	}

	public NFE fromUnitGeneratorVector(Vector<IntE> vector) {
		NFE result = one();
		List<NFE> unitGenerators = unitGroupGenerators();
		for (int i = 0; i < unitGenerators.size(); i++) {
			result = multiply(power(unitGenerators.get(i), vector.get(i + 1)), result);
		}
		return result;
	}

	public GenericPIDModule<IntE, Vector<IntE>> unitGroupAsIntModule() {
		FreeModule<IntE> module = new FreeModule<>(Integers.z(),
				field.realEmbeddings().size() + field.complexEmbeddings().size());
		return GenericPIDModule.fromSyzygies(module, Collections.singletonList(
				module.scalarMultiply(primitiveRootsOfUnity().lastKey(), module.getUnitVector(1)).asList()));
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
			result.add(row);
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

	private List<Vector<NFE>> simplifiedIntegerModuleGenerators(List<Vector<NFE>> generators) {
		if (generators.isEmpty()) {
			return Collections.emptyList();
		}
		FreeModule<NFE> free = new FreeModule<>(this, generators.get(0).dimension());
		List<Vector<IntE>> asIntVectors = new ArrayList<>();
		for (Vector<NFE> generator : generators) {
			for (NFE integral : getModuleGenerators()) {
				asIntVectors.add(asIntegerVector(free.scalarMultiply(integral, generator)));
			}
		}
		List<Vector<IntE>> integerResult = Integers.z().simplifySubModuleGenerators(Matrix.fromColumns(asIntVectors));
		List<Vector<NFE>> result = new ArrayList<>();
		for (Vector<IntE> v : integerResult) {
			result.add(fromIntegerVector(v));
		}
		return result;
	}

//	public SmallestIntegerSolutionPreparation prepareSmallestIntegerSolution(List<Vector<NFE>> generators) {
//		return prepareSmallestIntegerSolution(generators, Collections.emptyList());
//	}
//
//	public SmallestIntegerSolutionPreparation prepareSmallestIntegerSolution(List<Vector<NFE>> generators,
//			NumberFieldIdeal modulo) {
//		int dimension = generators.get(0).dimension();
//		FreeModule<NFE> free = new FreeModule<>(this, dimension);
//		List<Vector<NFE>> latticeGenerators = new ArrayList<>();
//		for (Vector<NFE> basisVector : free.getBasis()) {
//			for (NFE generator : modulo.generators()) {
//				latticeGenerators.add(free.scalarMultiply(generator, basisVector));
//			}
//		}
//		return prepareSmallestIntegerSolution(generators, latticeGenerators);
//	}
//
//
//	public SmallestIntegerSolutionPreparation prepareSmallestIntegerSolution(List<Vector<NFE>> generators,
//			List<Vector<NFE>> modulo) {
//		return Integers.z().prepareSmallestIntegerSolution(asIntegerVectorList(generators),
//				asIntegerVectorList(modulo));
//	}
//
//	public Vector<NFE> smallestKernelVector(SmallestIntegerSolutionPreparation preparation) {
//		return fromIntegerVector(Integers.z().smallestKernelVector(preparation));
//	}
//
//	public Vector<NFE> smallestKernelVector(List<Vector<NFE>> generators) {
//		return smallestKernelVector(prepareSmallestIntegerSolution(generators));
//	}
//
//	public Vector<NFE> smallestKernelVector(List<Vector<NFE>> generators, NumberFieldIdeal modulus) {
//		return smallestKernelVector(prepareSmallestIntegerSolution(generators, modulus));
//	}
//
//	public Vector<NFE> smallestKernelVector(List<Vector<NFE>> generators, List<Vector<NFE>> modulus) {
//		return smallestKernelVector(prepareSmallestIntegerSolution(generators, modulus));
//	}
//	
//	public Vector<NFE> smallestIntegerSolution(Vector<NFE> target, SmallestIntegerSolutionPreparation preparation) {
//		return fromIntegerVector(Integers.z().smallestIntegerSolution(asIntegerVector(target), preparation));
//	}
//
//	public Vector<NFE> smallestIntegerSolution(List<Vector<NFE>> generators, Vector<NFE> target) {
//		return smallestIntegerSolution(target, prepareSmallestIntegerSolution(generators));
//	}
//
//	public Vector<NFE> smallestIntegerSolution(List<Vector<NFE>> generators, Vector<NFE> target,
//			NumberFieldIdeal modulo) {
//		return smallestIntegerSolution(target, prepareSmallestIntegerSolution(generators, modulo));
//	}
//
//	public Vector<NFE> smallestIntegerSolution(List<Vector<NFE>> generators, Vector<NFE> target,
//			List<Vector<NFE>> modulo) {
//		return smallestIntegerSolution(target, prepareSmallestIntegerSolution(generators, modulo));
//	}

	public static class SmallestNumberFieldIntegerSolutionPreparation {
		private Matrix<IntE> generatorMatrix;
		private MatrixModule<IntE> matrixModule;
		private RealNumberFieldIntegerLattice kernelLattice;
		private FreeModule<NFE> solutionSpace;
		private MathMap<Vector<NFE>, Vector<NFE>> projectionMap;

		private SmallestNumberFieldIntegerSolutionPreparation(NumberFieldIntegers order, List<Vector<NFE>> generators,
				List<Vector<NFE>> modulus) {
			int accuracy = 128;
			this.projectionMap = new FunctionMathMap<>(
					(Vector<NFE> t) -> new Vector<>(t.asList().subList(0, generators.size())));
			this.solutionSpace = new FreeModule<>(order, generators.size());
			List<Vector<IntE>> actualGenerators = new ArrayList<>();
			actualGenerators.addAll(order.asIntegerMatrix(Matrix.fromColumns(generators)).asColumnList());
			for (Vector<NFE> mod : modulus) {
				actualGenerators.add(order.asIntegerVector(mod));
			}
			this.generatorMatrix = Matrix.fromColumns(actualGenerators);
			this.matrixModule = generatorMatrix.getModule(Integers.z());
			List<Vector<IntE>> kernelBasis = matrixModule.kernelBasis(generatorMatrix);
			if (kernelBasis.size() == 0) {
				return;
			}
			List<Vector<NFE>> projectedKernelBasis = new ArrayList<>();
			for (Vector<IntE> kernelVector : kernelBasis) {
				projectedKernelBasis.add(projectionMap.evaluate(order.fromIntegerVector(kernelVector)));
			}
			projectedKernelBasis = order.simplifiedIntegerModuleGenerators(projectedKernelBasis);
			if (projectedKernelBasis.size() == 0) {
				return;
			}
			Reals r = Reals.r(accuracy);
			this.kernelLattice = new RealNumberFieldIntegerLattice(r, generators.size(), order, projectedKernelBasis);
		}

		private Vector<NFE> smallestKernelVector() {
			return kernelLattice.getModuleGenerators().get(0);
		}

	}

	public SmallestNumberFieldIntegerSolutionPreparation prepareSmallestIntegerSolution(List<Vector<NFE>> generators) {
		return prepareSmallestIntegerSolution(generators, Collections.emptyList());
	}

	public SmallestNumberFieldIntegerSolutionPreparation prepareSmallestIntegerSolution(List<Vector<NFE>> generators,
			NumberFieldIdeal modulus) {
		int dimension = generators.get(0).dimension();
		FreeModule<NFE> free = new FreeModule<>(this, dimension);
		List<Vector<NFE>> latticeGenerators = new ArrayList<>();
		for (Vector<NFE> basisVector : free.getBasis()) {
			for (NFE generator : modulus.asSubModule.getBasis()) {
				latticeGenerators.add(free.scalarMultiply(generator, basisVector));
			}
		}
		return prepareSmallestIntegerSolution(generators, latticeGenerators);
	}

	public SmallestNumberFieldIntegerSolutionPreparation prepareSmallestIntegerSolution(List<Vector<NFE>> generators,
			List<Vector<NFE>> modulus) {
		return new SmallestNumberFieldIntegerSolutionPreparation(this, generators, modulus);
	}

	public SmallestNumberFieldIntegerSolutionPreparation prepareSmallestIntegerSolution(List<Vector<NFE>> generators,
			RealNumberFieldIntegerLattice modulus) {
		return prepareSmallestIntegerSolution(generators, modulus.getModuleGenerators());
	}

	public Vector<NFE> smallestKernelVector(SmallestNumberFieldIntegerSolutionPreparation preparation) {
		return preparation.smallestKernelVector();
	}

	public Vector<NFE> smallestKernelVector(List<Vector<NFE>> generators) {
		return smallestKernelVector(prepareSmallestIntegerSolution(generators));
	}

	public Vector<NFE> smallestKernelVector(List<Vector<NFE>> generators, NumberFieldIdeal modulus) {
		return smallestKernelVector(prepareSmallestIntegerSolution(generators, modulus));
	}

	public Vector<NFE> smallestKernelVector(List<Vector<NFE>> generators, List<Vector<NFE>> modulus) {
		return smallestKernelVector(prepareSmallestIntegerSolution(generators, modulus));
	}

	public Vector<NFE> smallestIntegerSolution(Vector<NFE> target,
			SmallestNumberFieldIntegerSolutionPreparation preparation) {
		Vector<NFE> solution = fromIntegerVector(
				preparation.matrixModule.solve(preparation.generatorMatrix, asIntegerVector(target)));
		if (preparation.kernelLattice == null) {
			return solution;
		}
		Vector<NFE> projectedSolution = preparation.projectionMap.evaluate(solution);
		Vector<NFE> closestKernelVector = preparation.kernelLattice.closestLatticePoint(projectedSolution);
		return preparation.solutionSpace.subtract(solution, closestKernelVector);
	}

	public Vector<NFE> smallestIntegerSolution(List<Vector<NFE>> generators, Vector<NFE> target) {
		return smallestIntegerSolution(target, prepareSmallestIntegerSolution(generators));
	}

	public Vector<NFE> smallestIntegerSolution(List<Vector<NFE>> generators, Vector<NFE> target,
			NumberFieldIdeal modulus) {
		return smallestIntegerSolution(target, prepareSmallestIntegerSolution(generators, modulus));
	}

	public Vector<NFE> smallestIntegerSolution(List<Vector<NFE>> generators, Vector<NFE> target,
			List<Vector<NFE>> modulus) {
		return smallestIntegerSolution(target, prepareSmallestIntegerSolution(generators, modulus));
	}

	@Override
	public List<Ideal<NFE>> maximalPrimeIdealChain(Ideal<NFE> start) {
		if (start.equals(getZeroIdeal())) {
			return maximalPrimeIdealChain(start,
					idealsOver(Integers.z().getIdeal(Collections.singletonList(Integers.z().getInteger(2)))).get(0));
		}
		return Collections.singletonList(primaryDecomposition(start).getRadicals().get(0));
	}

	@Override
	public List<Ideal<NFE>> maximalPrimeIdealChain(Ideal<NFE> start, Ideal<NFE> end) {
		if (!end.contains(start) || !end.isPrime()) {
			throw new ArithmeticException("Invalid preconditions!");
		}
		if (start.equals(getZeroIdeal())) {
			if (end.equals(getZeroIdeal())) {
				return Collections.singletonList(getZeroIdeal());
			}
			List<Ideal<NFE>> result = new ArrayList<>();
			result.add(start);
			result.add(end);
			return result;
		}
		return Collections.singletonList(end);
	}

	@Override
	public IdealResult<NFE, NumberFieldIdeal> getIdealWithTransforms(List<NFE> generators) {
		NumberFieldIdeal ideal = getIdeal(generators);
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
	public NumberFieldIdeal getIdeal(List<NFE> generators) {
		Integers z = Integers.z();
		IntE normGcd = z.zero();
		for (NFE generator : generators) {
			normGcd = z.gcd(normGcd, field.norm(generator).asInteger());
		}
		if (normGcd.equals(z.zero())) {
			List<NFE> expression = new ArrayList<>();
			for (int i = 0; i < generators.size(); i++) {
				expression.add(zero());
			}
			return getZeroIdeal();
		}
		FactorizationResult<IntE, IntE> factorization = z.uniqueFactorization(normGcd);
		return getIdealWithFactorization(generators, factorization);
	}

	@Override
	public NumberFieldIdeal getIdeal(NFE... generators) {
		return getIdeal(Arrays.asList(generators));
	}

	public Optional<NumberFieldIdeal> getIdealIfSmoothOver(List<NFE> generators, Set<IntE> smoothnessBase) {
		Integers z = Integers.z();
		IntE normGcd = z.zero();
		for (NFE generator : generators) {
			normGcd = z.gcd(normGcd, field.norm(generator).asInteger());
		}
		if (normGcd.equals(z.zero())) {
			return Optional.of(getZeroIdeal());
		}
		SortedMap<IntE, Integer> factors = new TreeMap<>();
		for (IntE prime : smoothnessBase) {
			int power = 0;
			while (z.isDivisible(normGcd, prime)) {
				normGcd = z.divide(normGcd, prime);
				power++;
			}
			if (power != 0) {
				factors.put(prime, power);
			}
		}
		if (!z.isUnit(normGcd)) {
			return Optional.empty();
		}
		return Optional.of(getIdealWithFactorization(generators, new FactorizationResult<>(normGcd, factors)));
	}

	private NumberFieldIdeal getIdealWithFactorization(List<NFE> generators,
			FactorizationResult<IntE, IntE> normFactorization) {
		if (normFactorization.primeFactors().isEmpty()) {
			return getUnitIdeal();
		}
		SortedMap<IntE, FactorizationResult<Ideal<NFE>, Ideal<NFE>>> factorizations = new TreeMap<>();
		SortedMap<Ideal<NFE>, Integer> allIdealFactors = new TreeMap<>();
		for (IntE prime : normFactorization.primeFactors()) {
			SortedMap<Ideal<NFE>, Integer> idealFactors = new TreeMap<>();
			List<NumberFieldIdeal> ideals = idealsOver(prime);
			for (NumberFieldIdeal ideal : ideals) {
				Value value = Value.INFINITY;
				for (NFE generator : generators) {
					if (generator.equals(zero())) {
						continue;
					}
					value = value.min(ideal.maximumPowerContains(generator));
				}
				if (!value.equals(Value.ZERO)) {
					idealFactors.put(ideal, value.value());
					allIdealFactors.put(ideal, value.value());
				}
			}
			if (normFactorization.primeFactors().size() == 1 && idealFactors.size() == 1
					&& idealFactors.get(idealFactors.firstKey()) == 1) {
				return (NumberFieldIdeal) idealFactors.firstKey();
			}
			factorizations.put(prime, new FactorizationResult<>(getUnitIdeal(), idealFactors));
		}
		FactorizationResult<Ideal<NFE>, Ideal<NFE>> factorization = new FactorizationResult<>(getUnitIdeal(),
				allIdealFactors);
		if (!idealsByFactorization.containsKey(factorization)) {
			if (factorization.isIrreducible()) {
				idealsByFactorization.put(factorization, (NumberFieldIdeal) factorization.firstPrimeFactor());
			} else {
				idealsByFactorization.put(factorization, new NumberFieldIdeal(factorizations));
			}
		}
		return idealsByFactorization.get(factorization);
	}

	@Override
	public NumberFieldIdeal multiply(Ideal<NFE> t1, Ideal<NFE> t2) {
		NumberFieldIdeal ideal1 = (NumberFieldIdeal) t1;
		NumberFieldIdeal ideal2 = (NumberFieldIdeal) t2;
		if (ideal1.zeroIdeal || ideal2.zeroIdeal) {
			return getZeroIdeal();
		}
		Set<Ideal<NFE>> primeFactors = new TreeSet<>();
		primeFactors.addAll(ideal1.idealFactorization.primeFactors());
		primeFactors.addAll(ideal2.idealFactorization.primeFactors());
		Map<NumberFieldIdeal, Integer> result = new TreeMap<>();
		for (Ideal<NFE> ideal : primeFactors) {
			result.put((NumberFieldIdeal) ideal,
					ideal1.idealFactorization.multiplicity(ideal) + ideal2.idealFactorization.multiplicity(ideal));
		}
		return fromFactorization(result);
	}

	@Override
	public NumberFieldIdeal intersect(Ideal<NFE> t1, Ideal<NFE> t2) {
		NumberFieldIdeal ideal1 = (NumberFieldIdeal) t1;
		NumberFieldIdeal ideal2 = (NumberFieldIdeal) t2;
		if (ideal1.zeroIdeal || ideal2.zeroIdeal) {
			return getZeroIdeal();
		}
		Set<Ideal<NFE>> primeFactors = new TreeSet<>();
		primeFactors.addAll(ideal1.idealFactorization.primeFactors());
		primeFactors.addAll(ideal2.idealFactorization.primeFactors());
		Map<NumberFieldIdeal, Integer> result = new TreeMap<>();
		for (Ideal<NFE> ideal : primeFactors) {
			result.put((NumberFieldIdeal) ideal, Math.max(ideal1.idealFactorization.multiplicity(ideal),
					ideal2.idealFactorization.multiplicity(ideal)));
		}
		return fromFactorization(result);
	}

	@Override
	public NumberFieldIdeal radical(Ideal<NFE> t) {
		NumberFieldIdeal ideal = (NumberFieldIdeal) t;
		if (ideal.zeroIdeal) {
			return getZeroIdeal();
		}
		Map<NumberFieldIdeal, Integer> result = new TreeMap<>();
		for (Ideal<NFE> idealFactor : ideal.idealFactorization.primeFactors()) {
			result.put((NumberFieldIdeal) idealFactor, 1);
		}
		return fromFactorization(result);
	}

	@Override
	public NumberFieldIdeal power(Ideal<NFE> t, int power) {
		Map<NumberFieldIdeal, Integer> factors = new TreeMap<>();
		FactorizationResult<Ideal<NFE>, Ideal<NFE>> originalFactors = idealFactorization(t);
		for (Ideal<NFE> prime : originalFactors.primeFactors()) {
			factors.put((NumberFieldIdeal) prime, power * originalFactors.multiplicity(prime));
		}
		return fromFactorization(factors);
	}

	NumberFieldIdeal fromFactorization(Map<NumberFieldIdeal, Integer> factors) {
		Map<NumberFieldIdeal, Integer> newFactors = new TreeMap<>();
		for (NumberFieldIdeal ideal : factors.keySet()) {
			int multiplicity = factors.get(ideal);
			if (multiplicity != 0) {
				newFactors.put(ideal, multiplicity);
			}
		}
		factors = newFactors;
		if (factors.size() == 0) {
			return getUnitIdeal();
		}
		if (factors.size() == 1 && factors.get(factors.keySet().iterator().next()) == 1) {
			return factors.keySet().iterator().next();
		}
		SortedMap<Ideal<NFE>, Integer> cast = new TreeMap<>();
		Map<IntE, SortedMap<Ideal<NFE>, Integer>> perPrimeFactorization = new TreeMap<>();
		for (Ideal<NFE> ideal : factors.keySet()) {
			NumberFieldIdeal primeIdeal = (NumberFieldIdeal) ideal;
			IntE prime = primeIdeal.intGenerator;
			if (!perPrimeFactorization.containsKey(prime)) {
				perPrimeFactorization.put(prime, new TreeMap<>());
			}
			perPrimeFactorization.get(prime).put(ideal, factors.get(ideal));
			cast.put(ideal, factors.get(ideal));
		}
		FactorizationResult<Ideal<NFE>, Ideal<NFE>> asFactorizationResult = new FactorizationResult<>(getUnitIdeal(),
				cast);
		if (!idealsByFactorization.containsKey(asFactorizationResult)) {
			Map<IntE, FactorizationResult<Ideal<NFE>, Ideal<NFE>>> perPrime = new TreeMap<>();
			for (IntE prime : perPrimeFactorization.keySet()) {
				perPrime.put(prime, new FactorizationResult<>(getUnitIdeal(), perPrimeFactorization.get(prime)));
			}
			idealsByFactorization.put(asFactorizationResult, new NumberFieldIdeal(perPrime));
		}
		return idealsByFactorization.get(asFactorizationResult);
	}

	@Override
	public FactorizationResult<Ideal<NFE>, Ideal<NFE>> idealFactorization(Ideal<NFE> t) {
		NumberFieldIdeal ideal = (NumberFieldIdeal) t;
		return ideal.idealFactorization;
	}

	@Override
	public DedekindRing<NFE, NFE, FFE> asDedekindRing() {
		return this;
	}

	@Override
	public PrimaryDecompositionResult<NFE, NumberFieldIdeal> primaryDecomposition(Ideal<NFE> t) {
		NumberFieldIdeal ideal = (NumberFieldIdeal) t;
		List<NumberFieldIdeal> primaries = new ArrayList<>();
		List<NumberFieldIdeal> radicals = new ArrayList<>();
		for (Ideal<NFE> prime : ideal.idealFactorization.primeFactors()) {
			NumberFieldIdeal primeIdeal = (NumberFieldIdeal) prime;
			primaries.add(power(prime, ideal.idealFactorization.multiplicity(prime)));
			radicals.add(primeIdeal);
		}
		return new PrimaryDecompositionResult<>(primaries, radicals);
	}

	@Override
	public ModuloMaximalIdealResult<NFE, FFE, NumberFieldIntegers, NumberFieldIdeal, FiniteField> moduloMaximalIdeal(
			Ideal<NFE> ideal) {
		return new ModuloMaximalIdealResult<>(this, (NumberFieldIdeal) ideal, reduction(ideal), new MathMap<>() {
			@Override
			public FFE evaluate(NFE t) {
				return reduce(t, ideal);
			}
		}, new MathMap<>() {
			@Override
			public NFE evaluate(FFE t) {
				return lift(t, ideal);
			}
		});
	}

	@Override
	public ModuloIdealResult<NFE, ?> moduloIdeal(Ideal<NFE> ideal) {
		if (ideal.isMaximal()) {
			ModuloMaximalIdealResult<NFE, FFE, NumberFieldIntegers, NumberFieldIdeal, FiniteField> result = moduloMaximalIdeal(
					ideal);
			return new ModuloIdealResult<>(this, ideal, result.getField(), result.getReduction(), result.getLift());
		}
		if (ideal.isPrime()) {
			return new ModuloIdealResult<>(this, ideal, this, new Identity<>(), new Identity<>());
		}
		ModuloNumberFieldIdeal modulo = ((NumberFieldIdeal) ideal).modOut();// new ModuloNumberFieldIdeal(this,
																			// (NumberFieldIdeal) ideal);
		return new ModuloIdealResult<>(this, ideal, modulo, new MathMap<>() {
			@Override
			public ModNFE evaluate(NFE t) {
				return modulo.reduce(t);
			}

		}, new MathMap<>() {
			@Override
			public NFE evaluate(ModNFE t) {
				return modulo.lift(t);
			}
		});
	}

	@Override
	public NumberFieldIdeal getUnitIdeal() {
		if (unitIdeal == null) {
			this.unitIdeal = new NumberFieldIdeal(false);
			this.idealsByFactorization.put(unitIdeal.idealFactorization, unitIdeal);
		}
		return unitIdeal;
	}

	@Override
	public NumberFieldIdeal getZeroIdeal() {
		if (zeroIdeal == null) {
			this.zeroIdeal = new NumberFieldIdeal(true);
		}
		return zeroIdeal;
	}

	@Override
	public NumberFieldIdeal getNilRadical() {
		return getZeroIdeal();
	}

	public List<NumberFieldIdeal> idealsOver(int prime) {
		return idealsOver(BigInteger.valueOf(prime));
	}

	public List<NumberFieldIdeal> idealsOver(BigInteger prime) {
		return idealsOver(Integers.z().getInteger(prime));
	}

	public List<NumberFieldIdeal> idealsOver(Ideal<IntE> integerIdeal) {
		return idealsOver(integerIdeal.generators().get(0));
	}

	public List<NumberFieldIdeal> idealsOver(IntE prime) {
		Integers z = Integers.z();
		if (!z.isPrime(prime)) {
			throw new ArithmeticException("Prime number required!");
		}
		if (idealsOverPrime.containsKey(prime)) {
			return idealsOverPrime.get(prime);
		}
//		DiscreteValuationRing<Fraction, PFE> zp = Integers.z().localize(prime);
//		PrimeField fp = PrimeField.getPrimeField(prime);
//		Rationals q = Rationals.q();
//		UnivariatePolynomial<Fraction> zpMinimalPolynomial = q.getUnivariatePolynomialRing()
//				.getEmbedding(minimalPolynomial(), new MathMap<>() {
//					@Override
//					public Fraction evaluate(IntE t) {
//						return q.getInteger(t);
//					}
//				});
//		TheMontesResult<Fraction, PFE, PFE, FFE, FiniteField> theMontes = zp.theMontesAlgorithm(zpMinimalPolynomial,
//				fp.getExtension(fp.getUnivariatePolynomialRing().getVar()));
//		List<NumberFieldIdeal> ideals = new ArrayList<>();
//		List<List<UnivariatePolynomial<Fraction>>> integralBasis = zp.integralBasis(zpMinimalPolynomial, theMontes,
//				true);
//		for (int i = 0; i < theMontes.getTypes().size(); i++) {
//			ideals.add(new NumberFieldIdeal(zp, theMontes.getTypes(), i, integralBasis));
//		}
		List<TwoGeneratorIdeal<IntE, Fraction, PFE, PFE, FFE, FiniteField>> twoGeneratorIdeals = z
				.idealsOverGenerators(minimalPolynomial(), prime, PrimeField.getPrimeField(prime)
						.getExtension(PrimeField.getPrimeField(prime).getUnivariatePolynomialRing().getVar()));
		List<NumberFieldIdeal> ideals = new ArrayList<>();
		for (TwoGeneratorIdeal<IntE, Fraction, PFE, PFE, FFE, FiniteField> ideal : twoGeneratorIdeals) {
			ideals.add(new NumberFieldIdeal(ideal));
		}
		idealsOverPrime.put(prime, ideals);
		return ideals;
	}

	public List<NumberFieldIdeal> idealsOver(Ideal<NFE> lowerIdeal,
			FieldEmbedding<Fraction, NFE, NumberField> fieldEmbedding) {
		if (!fieldEmbedding.getField().equals(field)) {
			throw new ArithmeticException("incorrect field embedding");
		}
		NumberFieldIdeal lowerPrimeIdeal = (NumberFieldIdeal) lowerIdeal;
		NFE embeddedLowerUniformizer = fieldEmbedding.getEmbedding(lowerPrimeIdeal.uniformizer);
		Integers z = Integers.z();
		List<NumberFieldIdeal> overRationals = idealsOver(z.getIdeal(lowerPrimeIdeal.intGenerator));
		List<NumberFieldIdeal> result = new ArrayList<>();
		for (NumberFieldIdeal ideal : overRationals) {
			if (ideal.contains(embeddedLowerUniformizer)) {
				result.add(ideal);
			}
		}
		return result;
	}

	public IntegerIdeal intersectToIntegers(Ideal<NFE> ideal) {
		return Integers.z().getIdeal(Collections.singletonList(((NumberFieldIdeal) ideal).intGenerator));
	}

	public NumberFieldIdeal intersectToLowerField(Ideal<NFE> ideal,
			FieldEmbedding<Fraction, NFE, NumberField> fieldEmbedding) {
		NumberFieldIdeal numberFieldIdeal = (NumberFieldIdeal) ideal;
		NumberFieldIntegers lowerOrder = fieldEmbedding.getEmbeddedField().maximalOrder();
		List<NFE> lowerIntegralBasis = new ArrayList<>();
		for (NFE t : lowerOrder.getModuleGenerators()) {
			lowerIntegralBasis.add(fieldEmbedding.getEmbedding(t));
		}
		FreeSubModule<IntE, NFE> lowerIntegers = new FreeSubModule<>(this, lowerIntegralBasis);
		FreeSubModule<IntE, NFE> intersected = lowerIntegers.intersection(numberFieldIdeal.asSubModule);
		List<NFE> lowerGenerators = new ArrayList<>();
		for (NFE intersect : intersected.getBasis()) {
			Vector<NFE> asVector = fieldEmbedding.asVector(intersect);
			for (int i = 1; i < asVector.dimension(); i++) {
				if (!asVector.get(i + 1).equals(lowerOrder.zero())) {
					throw new ArithmeticException("Could not intersect");
				}
			}
			lowerGenerators.add(asVector.get(1));
		}
		return lowerOrder.getIdeal(lowerGenerators);
	}

	public NumberFieldIdeal extend(Ideal<IntE> ideal) {
		return getIdeal(Collections.singletonList(getEmbedding(ideal.generators().get(0))));
	}

	public NumberFieldIdeal extend(Ideal<NFE> ideal, FieldEmbedding<Fraction, NFE, NumberField> fieldEmbedding) {
		if (!fieldEmbedding.getField().equals(field)) {
			throw new ArithmeticException("Field Embedding incorrect");
		}
		return (NumberFieldIdeal) getIdealEmbedding(ideal, fieldEmbedding.getEmbeddingMap());
	}

	public IntegerIdeal norm(Ideal<NFE> ideal) {
		NumberFieldIdeal numberFieldIdeal = (NumberFieldIdeal) ideal;
		return Integers.z().getIdeal(Collections.singletonList(numberFieldIdeal.norm()));
	}

	public NumberFieldIdeal normOver(Ideal<NFE> ideal, FieldEmbedding<Fraction, NFE, NumberField> fieldEmbedding) {
		if (!fieldEmbedding.getField().equals(field)) {
			throw new ArithmeticException("Field Embedding incorrect");
		}
		NumberFieldIntegers lowerOrder = fieldEmbedding.getEmbeddedField().maximalOrder();
		if (ideal.equals(getZeroIdeal())) {
			return lowerOrder.getZeroIdeal();
		}
		if (ideal.equals(getUnitIdeal())) {
			return lowerOrder.getUnitIdeal();
		}
		NumberFieldIdeal numberFieldIdeal = (NumberFieldIdeal) ideal;
		if (numberFieldIdeal.isMaximal()) {
			NumberFieldIdeal intersected = intersectToLowerField(numberFieldIdeal, fieldEmbedding);
			int degree = numberFieldIdeal.type().residueDegree() / intersected.type().residueDegree();
			return lowerOrder.power(intersected, degree);
		}
		NumberFieldIdeal result = lowerOrder.getUnitIdeal();
		for (Ideal<NFE> primeIdeal : numberFieldIdeal.idealFactorization.primeFactors()) {
			result = lowerOrder.multiply(result, lowerOrder.power(normOver(primeIdeal, fieldEmbedding),
					numberFieldIdeal.idealFactorization.multiplicity(primeIdeal)));
		}
		return result;
	}

	private UnivariatePolynomial<IntE> minimalPolynomial() {
		if (this.minimalPolynomial == null) {
			UnivariatePolynomial<Fraction> minimalPolynomial = field.minimalPolynomial();
			Integers z = Integers.z();
			Rationals q = Rationals.q();
			UnivariatePolynomialRing<Fraction> rationalPolynomials = q.getUnivariatePolynomialRing();
			if (!minimalPolynomial.leadingCoefficient().equals(q.one())) {
				throw new ArithmeticException("Minimal polynomial should be normalized");
			}
			IntE denom = z.one();
			for (int i = 0; i < minimalPolynomial.degree(); i++) {
				denom = z.lcm(denom, minimalPolynomial.univariateCoefficient(i).getDenominator());
			}
			minimalPolynomial = rationalPolynomials
					.normalize(rationalPolynomials.substitute(minimalPolynomial, Collections.singletonList(
							rationalPolynomials.divideScalar(rationalPolynomials.getVar(), q.getInteger(denom)))));
			this.minimalPolynomialDenominator = denom;
			this.minimalPolynomial = z.getUnivariatePolynomialRing().getEmbedding(minimalPolynomial, new MathMap<>() {
				@Override
				public IntE evaluate(Fraction t) {
					return t.getNumerator();
				}
			});
		}
		return minimalPolynomial;
	}

	public NFE different(NFE t) {
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> polynomials = q.getUnivariatePolynomialRing();
		UnivariatePolynomial<Fraction> minimalPolynomial = field.minimalPolynomial(t);
		if (minimalPolynomial.degree() != field.degree()) {
			return zero();
		}
		UnivariatePolynomial<Fraction> derivative = polynomials.derivative(minimalPolynomial);
		UnivariatePolynomialRing<NFE> numberFieldPolynomials = field.getUnivariatePolynomialRing();
		UnivariatePolynomial<NFE> embeddedDerivative = numberFieldPolynomials.getEmbedding(derivative,
				field.getEmbeddingMap());
		return numberFieldPolynomials.evaluate(embeddedDerivative, t);
	}

	public NFE differentOver(NFE t, FieldEmbedding<Fraction, NFE, NumberField> fieldEmbedding) {
		UnivariatePolynomialRing<NFE> polynomials = fieldEmbedding.getEmbeddedField().getUnivariatePolynomialRing();
		UnivariatePolynomial<NFE> minimalPolynomial = fieldEmbedding.minimalPolynomial(t);
		if (minimalPolynomial.degree() != fieldEmbedding.relativeDegree()) {
			return zero();
		}
		UnivariatePolynomial<NFE> derivative = polynomials.derivative(minimalPolynomial);
		UnivariatePolynomialRing<NFE> numberFieldPolynomials = fieldEmbedding.getField().getUnivariatePolynomialRing();
		UnivariatePolynomial<NFE> embeddedDerivative = numberFieldPolynomials.getEmbedding(derivative,
				fieldEmbedding.getEmbeddingMap());
		return numberFieldPolynomials.evaluate(embeddedDerivative, t);
	}

	public NumberFieldIdeal different() {
		List<NFE> differentGenerator = new ArrayList<>();
		for (IntE prime : potentiallyRamifiedPrimes()) {
			LocalizedFractions localized = Rationals.q().withValuation(prime.getValue());
			TheMontesResult<Fraction, PFE, PFE, FFE, FiniteField> theMontes = localized.ringOfIntegers()
					.theMontesAlgorithm(field.minimalPolynomial(), localized.residueField()
							.getExtension(localized.residueField().getUnivariatePolynomialRing().getVar()));
			for (OkutsuType<Fraction, PFE, PFE, FFE, FiniteField> type : theMontes.getTypes()) {
				for (int j = 0; j < type.representative().degree(); j++) {
					differentGenerator.add(different(field.fromPolynomial(type.divisorPolynomial(j))));
				}
			}
		}
		for (NFE generator : getModuleGenerators()) {
			for (int i = 1; i < field.degree(); i++) {
				differentGenerator.add(different(power(generator, i)));
			}
		}
		return getIdeal(differentGenerator);
	}

	public NumberFieldIdeal differentOver(FieldEmbedding<Fraction, NFE, NumberField> fieldEmbedding) {
		if (!field.equals(fieldEmbedding.getField())) {
			throw new ArithmeticException("needs to be called on upper order");
		}
		List<NFE> differentGenerator = new ArrayList<>();
		NumberFieldIntegers order = fieldEmbedding.getEmbeddedField().maximalOrder();
		for (IntE prime : order.potentiallyRamifiedPrimes()) {
			for (NumberFieldIdeal primeIdeal : order.idealsOver(prime)) {
				LocalizedNumberField localized = order.localizeAndQuotient(primeIdeal);
				TheMontesResult<NFE, FFE, PFE, FFE, FiniteField> theMontes = localized.ringOfIntegers()
						.theMontesAlgorithm(fieldEmbedding.minimalPolynomial(), localized.residueField()
								.getExtension(localized.residueField().getUnivariatePolynomialRing().getVar()));
				List<List<UnivariatePolynomial<NFE>>> integralBasis = localized.ringOfIntegers()
						.integralBasis(fieldEmbedding.minimalPolynomial(), theMontes, true);
				for (List<UnivariatePolynomial<NFE>> integralList : integralBasis) {
					for (UnivariatePolynomial<NFE> t : integralList) {
						differentGenerator.add(differentOver(fieldEmbedding.fromPolynomial(t), fieldEmbedding));
					}
				}
			}
		}
		for (NFE generator : getModuleGenerators()) {
			for (int i = 1; i < field.degree(); i++) {
				differentGenerator.add(differentOver(power(generator, i), fieldEmbedding));
			}
		}
		return getIdeal(differentGenerator);
	}

	public IntE discriminant() {
		return discriminant(getModuleGenerators());
	}

	IntE discriminant(List<NFE> generators) {
		if (!field.isBasis(generators)) {
			throw new ArithmeticException("Discriminant needs a basis");
		}
		List<List<Fraction>> traceForm = new ArrayList<>();
		for (int i = 0; i < generators.size(); i++) {
			if (!field.isInteger(generators.get(i))) {
				throw new ArithmeticException("Not an integer element");
			}
			List<Fraction> traceFormRow = new ArrayList<>();
			for (int j = 0; j < generators.size(); j++) {
				traceFormRow.add(field.trace(field.multiply(generators.get(i), generators.get(j))));
			}
			traceForm.add(traceFormRow);
		}
		Matrix<Fraction> trace = new Matrix<>(traceForm);
		return field.matrixAlgebra().determinant(trace).asInteger();
	}

	public NumberFieldIdeal discriminantOver(FieldEmbedding<Fraction, NFE, NumberField> fieldEmbedding) {
		return differentOver(fieldEmbedding).normOver(fieldEmbedding);
	}

	private Set<IntE> potentiallyRamifiedPrimes() {
		Integers z = Integers.z();
		return z.uniqueFactorization(z.getUnivariatePolynomialRing().discriminant(minimalPolynomial())).primeFactors();
	}

	public IntE norm(NFE t) {
		return field.norm(t).asInteger();
	}

	public IntE trace(NFE t) {
		return field.trace(t).asInteger();
	}

	@Override
	public Integers getRing() {
		return Integers.z();
	}

	@Override
	public boolean isFree() {
		return true;
	}

	@Override
	public IntegerIdeal annihilator() {
		return Integers.z().getZeroIdeal();
	}

	@Override
	public List<Vector<IntE>> getSyzygies() {
		return Collections.emptyList();
	}

	@Override
	public boolean isLinearIndependent(List<NFE> s) {
		List<Vector<IntE>> asVectors = new ArrayList<>();
		for (NFE t : s) {
			asVectors.add(asVector(t));
		}
		return asFreeModule.isLinearIndependent(asVectors);
	}

	@Override
	public boolean isGeneratingModule(List<NFE> s) {
		List<Vector<IntE>> asVectors = new ArrayList<>();
		for (NFE t : s) {
			asVectors.add(asVector(t));
		}
		return asFreeModule.isGeneratingModule(asVectors);
	}

	@Override
	public List<Vector<IntE>> nonTrivialCombinations(List<NFE> s) {
		List<Vector<IntE>> asVectors = new ArrayList<>();
		for (NFE t : s) {
			asVectors.add(asVector(t));
		}
		return asFreeModule.nonTrivialCombinations(asVectors);
	}

	public boolean isElement(NFE t) {
		Rationals q = Rationals.q();
		if (toIntegralBasis == null) {
			UnivariatePolynomial<Fraction> minPoly = field.minimalPolynomial(t);
			if (!minPoly.leadingCoefficient().equals(q.one())) {
				minPoly = q.getUnivariatePolynomialRing().normalize(minPoly);
			}
			for (int i = 0; i < minPoly.degree(); i++) {
				if (!q.isInteger(minPoly.univariateCoefficient(i))) {
					return false;
				}
			}
			return true;
		}
		Vector<Fraction> integerVector = field.matrixAlgebra().multiply(toIntegralBasis, field.asVector(t));
		for (Fraction c : integerVector.asList()) {
			if (!q.isInteger(c)) {
				return false;
			}
		}
		return true;
	}

	@Override
	public List<NFE> getModuleGenerators() {
		if (integralBasis == null) {
			Rationals q = Rationals.q();
			Integers z = Integers.z();
			if (field.degree() == 1) {
				this.integralBasis = Collections.singletonList(one());
				this.fromIntegralBasis = field.matrixAlgebra().one();
				this.toIntegralBasis = field.matrixAlgebra().one();
				return integralBasis;
			}
			UnivariatePolynomialRing<Fraction> polynomials = q.getUnivariatePolynomialRing();
			List<IntE> primes = new ArrayList<>();
			primes.addAll(potentiallyRamifiedPrimes());
			Map<IntE, List<UnivariatePolynomial<Fraction>>> integralBasisPerPrime = new TreeMap<>();
			Map<IntE, List<Integer>> valuesPerPrime = new TreeMap<>();
			BigInteger allPrimes = BigInteger.ONE;
			UnivariatePolynomial<Fraction> zpMinimalPolynomial = q.getUnivariatePolynomialRing()
					.getEmbedding(minimalPolynomial(), new MathMap<>() {
						@Override
						public Fraction evaluate(IntE t) {
							return q.getInteger(t);
						}
					});
			for (IntE prime : primes) {
				allPrimes = allPrimes.multiply(prime.getValue());
				PrimeField fp = PrimeField.getPrimeField(prime.getValue());
				DiscreteValuationRing<Fraction, PFE> zp = z.localize(prime.getValue());
				TheMontesResult<Fraction, PFE, PFE, FFE, FiniteField> theMontes = zp.theMontesAlgorithm(
						zpMinimalPolynomial, fp.getExtension(fp.getUnivariatePolynomialRing().getVar()));
				List<UnivariatePolynomial<Fraction>> integral = zp.triagonalizeIntegralBasis(zpMinimalPolynomial,
						zp.integralBasis(zpMinimalPolynomial, theMontes, false));
				List<Integer> valueList = new ArrayList<>();
				for (int i = 0; i < minimalPolynomial().degree(); i++) {
					valueList.add(-zp.localField().valuation(integral.get(i).leadingCoefficient()).value());
				}
				integralBasisPerPrime.put(prime, integral);
				valuesPerPrime.put(prime, valueList);
			}
			this.integralBasis = new ArrayList<>();
			for (int i = 0; i < minimalPolynomial().degree(); i++) {
				UnivariatePolynomial<Fraction> b = polynomials.zero();
				IntE primeProduct = z.one();
				IntE primeProductDivisor = z.one();
				for (IntE prime : primes) {
					IntE primePower = z.power(prime, valuesPerPrime.get(prime).get(i));
					primeProductDivisor = z.multiply(primePower, primeProductDivisor);
					primeProduct = z.multiply(prime, primePower, primeProduct);
				}
				for (IntE prime : primes) {
					BigInteger primePower = z.power(prime, valuesPerPrime.get(prime).get(i) + 1).getValue();
					IntE allOthers = z.divideChecked(primeProduct, z.getInteger(prime));
					IntE multiplier = z.multiply(primeProduct.getValue().divide(primePower).modInverse(primePower),
							allOthers);
					b = polynomials.add(b, polynomials.multiply(multiplier, integralBasisPerPrime.get(prime).get(i)));
				}
				BigInteger product = primeProduct.getValue();
				b = polynomials.getEmbedding(b, new MathMap<>() {
					@Override
					public Fraction evaluate(Fraction t) {
						BigInteger rounded = t.getNumerator().getValue().divide(t.getDenominator().getValue());
						BigInteger quotient = rounded.divide(product);
						return q.subtract(t, q.getInteger(quotient.multiply(product)));
					}
				});
				b = polynomials.divideScalar(b, q.getEmbedding(primeProductDivisor));
				if (b.degree() != i) {
					throw new ArithmeticException("expected integral basis to be triagonal!");
				}

				NFE element = field.fromPolynomial(polynomials.substitute(b, Collections
						.singletonList(polynomials.getEmbedding(q.getInteger(minimalPolynomialDenominator), 1))));
//				if (!isElement(element)) {
//					throw new ArithmeticException("Integral basis algorithm wrong, " + element + " is not integer");
//				}
				integralBasis.add(element);
			}
			this.integralBasis = sublatticeReduction(integralBasis);
			List<Vector<Fraction>> asVectors = new ArrayList<>();
			for (NFE b : integralBasis) {
				asVectors.add(field.asVector(b));
			}
			this.fromIntegralBasis = Matrix.fromColumns(asVectors);
			this.toIntegralBasis = field.matrixAlgebra().inverse(fromIntegralBasis);
		}
		return integralBasis;
	}

	@Override
	public Vector<IntE> asVector(NFE s) {
		Vector<Fraction> asFieldVector = field.asVector(s);
		Vector<Fraction> asIntVector = field.matrixAlgebra().multiply(toIntegralBasisBaseChange(), asFieldVector);
		List<IntE> ints = new ArrayList<>();
		for (Fraction c : asIntVector.asList()) {
			ints.add(c.asInteger());
		}
		return new Vector<>(ints);
	}

	public NFE roundToInteger(NFE t) {
		Vector<Fraction> inIntegralBasis = field.matrixAlgebra().multiply(toIntegralBasisBaseChange(),
				field.asVector(t));
		List<IntE> rounded = new ArrayList<>();
		for (Fraction coefficient : inIntegralBasis.asList()) {
			rounded.add(coefficient.round());
		}
		return fromVector(new Vector<>(rounded));
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

	public Matrix<Fraction> toIntegralBasisBaseChange() {
		getModuleGenerators();
		return toIntegralBasis;
	}

	@Override
	public Value valuation(NFE t, Ideal<NFE> maximalIdeal) {
		return maximalIdeal.maximumPowerContains(t);
	}

	public NFE getNumerator(NFE t) {
		if (isElement(t)) {
			return t;
		}
		return field.multiply(t, getDenominator(t));
	}

	public NFE getDenominator(NFE t) {
		return getEmbedding(getIntegerDenominator(t));
	}

	public IntE getIntegerDenominator(NFE t) {
		Integers z = Integers.z();
		if (isElement(t)) {
			return z.one();
		}
		Vector<Fraction> asVector = field.asVector(t);
		Matrix<Fraction> baseChange = toIntegralBasisBaseChange();
		Vector<Fraction> asIntegralVector = field.asVectorSpace().matrixAlgebra().multiply(baseChange, asVector);
		IntE denominator = z.one();
		for (Fraction c : asIntegralVector.asList()) {
			denominator = z.lcm(denominator, c.getDenominator());
		}
		return denominator;
	}

	@Override
	public NFE asInteger(NFE t) {
		if (isElement(t)) {
			return t;
		}
		throw new ArithmeticException("Not an integer!");
	}

	@Override
	public FieldOfFractionsResult<NFE, NFE> fieldOfFractions() {
		return new FieldOfFractionsResult<>(this, field, new Identity<>(), new MathMap<>() {

			@Override
			public NFE evaluate(NFE t) {
				return getNumerator(t);
			}
		}, new MathMap<>() {

			@Override
			public NFE evaluate(NFE t) {
				return getDenominator(t);
			}
		}, new MathMap<>() {

			@Override
			public NFE evaluate(NFE t) {
				return asInteger(t);
			}
		});
	}

	// TODO: Fix
	@Override
	public LocalizeResult<NFE, NFE, NFE, FFE> localizeAtIdeal(Ideal<NFE> primeIdeal) {
		localize(primeIdeal);
		LocalizedNumberField localized = localizations.get(primeIdeal);
		return new LocalizeResult<>(this, primeIdeal, localized.ringOfIntegers(), new Identity<>(), new MathMap<>() {
			@Override
			public NFE evaluate(NFE t) {
				Value v = localized.valuation(t);
				if (v.compareTo(Value.ZERO) >= 0) {
					return t;
				}
				return field.multiply(t, power(localized.uniformizer(), -v.value()));
			}
		}, new MathMap<>() {

			@Override
			public NFE evaluate(NFE t) {
				Value v = localized.valuation(t);
				if (v.compareTo(Value.ZERO) >= 0) {
					return one();
				}
				return power(localized.uniformizer(), -v.value());
			}
		}, new MathMap<>() {

			@Override
			public NFE evaluate(NFE t) {
				if (!isElement(t)) {
					throw new ArithmeticException("Not an integer!");
				}
				return t;
			}
		});
	}

	@Override
	public boolean isInteger(NFE t) {
		return isElement(t);
	}

	@Override
	public NumberField quotientField() {
		return numberField();
	}

	public NumberField numberField() {
		return field;
	}

	@Override
	public FiniteRealVectorSpace getVectorSpace() {
		return field.minkowskiEmbeddingSpace();
	}

	@Override
	public Vector<Real> embedding(NFE t) {
		return field.minkowskiEmbedding(t);
	}

	List<NFE> sublatticeReduction(List<NFE> sublatticeBasis) {
		return getVectorSpace().latticeReduction(sublatticeBasis, this);
	}

	@Override
	public int rank() {
		return field.degree();
	}

	public NumberFieldOrder asOrder() {
		return field.getOrder(getModuleGenerators());
	}

	FractionalIdeal reducedRepresentative(Ideal<NFE> b) {
		IdealGroup idealGroup = field.idealGroup();
		FractionalIdeal ideal = idealGroup.getEmbedding(b);
		NFE minimal = ((NumberFieldIdeal) b).shortestElement();
		FractionalIdeal result = idealGroup.operate(ideal, idealGroup.getPrincipalIdeal(field.inverse(minimal)));
		return result;
	}

	List<NumberFieldIdeal> idealEquivalenceClass(Ideal<NFE> b) {
		boolean calculateUnits = units == null && field.realEmbeddings().size() + field.complexEmbeddings().size() >= 2
				&& b.equals(getUnitIdeal());
		int k = 0;
		List<FractionalIdeal> idealList = new ArrayList<>();
		idealList.add(reducedRepresentative(b));
		List<NFE> multipliers = null;
		if (calculateUnits) {
			multipliers = new ArrayList<>();
			multipliers.add(one());
			units = new TreeSet<>();
		}
		IdealGroup idealGroup = field.idealGroup();
		while (k < idealList.size()) {
			for (NFE t : neighbors(idealList.get(k))) {
				FractionalIdeal nextIdeal = idealGroup.operate(idealList.get(k),
						idealGroup.getPrincipalIdeal(field.inverse(t)));
				boolean found = false;
				for (int j = 0; j < idealList.size(); j++) {
					if (nextIdeal.equals(idealList.get(j))) {
						if (calculateUnits) {
							NFE newUnit = field.multiply(t, field.divide(multipliers.get(k), multipliers.get(j)));
							if (!isRootOfUnity(newUnit)) {
								units.add(newUnit);
							}
						}
						found = true;
						break;
					}
				}
				if (!found) {
					idealList.add(nextIdeal);
					if (calculateUnits) {
						multipliers.add(field.multiply(t, multipliers.get(k)));
					}
				}
			}
			k++;
		}
		List<NumberFieldIdeal> ideals = new ArrayList<>();
		for (FractionalIdeal ideal : idealList) {
			ideals.add(ideal.clearDenominator().getFirst());
		}
		if (calculateUnits) {
			this.unitIdealMultipliers = multipliers;
			List<NFE> unitList = new ArrayList<>();
			unitList.addAll(units);
			units.clear();
			List<Vector<Real>> logarithmList = new ArrayList<>();
			for (NFE unit : unitList) {
				logarithmList.add(field.logRepresentation(unit));
			}
			unitLattice = new RealLattice(field.logRepresentationSpace(), logarithmList);
			if (unitLattice.rank() != field.realEmbeddings().size() + field.complexEmbeddings().size() - 1) {
				throw new ArithmeticException("Unit computation went wrong!");
			}
			List<Vector<IntE>> asIntVectors = new ArrayList<>();
			for (Vector<Real> logVector : logarithmList) {
				asIntVectors.add(unitLattice.asVector(logVector));
			}
			Matrix<IntE> asMatrix = Matrix.fromColumns(asIntVectors);
			MatrixModule<IntE> module = asMatrix.getModule(Integers.z());
			FreeModule<IntE> free = new FreeModule<>(Integers.z(), unitLattice.rank());
			for (int i = 0; i < unitLattice.rank(); i++) {
				Vector<IntE> exponents = module.solve(asMatrix, free.getUnitVector(i + 1));
				NFE unit = one();
				for (int j = 0; j < unitList.size(); j++) {
					unit = multiply(power(unitList.get(j), exponents.get(j + 1)), unit);
				}
				units.add(unit);
			}
		}
		return ideals;
	}

	private boolean isRootOfUnity(NFE unit) {
		Reals r = field.minkowskiEmbeddingSpace().getValueField();
		for (Real value : values(unit).asList()) {
			if (r.abs(r.subtract(value, r.one())).compareTo(r.getPowerOfTwo(-r.precision() / 2)) > 0) {
				return false;
			}
		}
		return true;
	}

	private Set<NFE> neighbors(FractionalIdeal ideal) {
		int l = 0;
		Reals r = field.minkowskiEmbeddingSpace().getValueField();
		List<Set<NFE>> minimalSets = new ArrayList<>();
		Set<NFE> startingSet = new TreeSet<>();
		List<Real> border = new ArrayList<>();
		List<Dual<Real, Vector<Real>>> constraints = new ArrayList<>();
		int realEmbeddings = field.realEmbeddings().size();
		int complexEmbeddings = field.complexEmbeddings().size();
		int numComplexConstraints = 1 << (field.degree() + 1);
		constraints = circle(numComplexConstraints, -1);
		for (int i = 0; i < realEmbeddings + complexEmbeddings; i++) {
			if (i < realEmbeddings) {
				border.add(r.one());
				border.add(r.one());
			} else {
				for (int j = 0; j < numComplexConstraints; j++) {
					border.add(r.one());
				}
			}
		}
		Real epsilon = r.getPowerOfTwo(-r.precision() / 2);
		Real boundary = r.add(r.one(), epsilon);
//		System.err.println(ideal);
		candidateLoop: for (NFE candidate : field.minkowskiEmbeddingSpace()
				.latticePointsInPolytope(new Polytope<>(field.minkowskiEmbeddingSpace(), constraints, border), ideal)) {
			if (candidate.equals(zero())) {
				continue;
			}
			Vector<Real> values = values(candidate);
			for (int i = 0; i < values.dimension(); i++) {
				if (values.get(i + 1).compareTo(boundary) > 0) {
					continue candidateLoop;
				}
			}
			startingSet.add(candidate);
		}
//		System.err.println("Ideal: " + ideal);
//		System.err.println("Starting Set: " + startingSet);
		minimalSets.add(startingSet);
		int embeddings = field.logRepresentationSpace().dimension();
		while (l < minimalSets.size()) {
//			System.out.println("l: " + l + "\nMinimal Sets: " + minimalSets);
			for (int i = 0; i < embeddings; i++) {
				Set<NFE> set = minimalSets.get(l);
				Set<NFE> expansion = expansion(set, i + 1, ideal);
//				System.out.println("Expansion " + (i + 1) + ": " + expansion);
				if (!minimalSets.contains(expansion)) {
//					System.err.println("Added Expansion: " + expansion);
					minimalSets.add(expansion);
				}
				if (value(set, i + 1).compareTo(r.one()) > 0) {
					Set<NFE> compression = compression(set, i + 1);
//					System.out.println("Compression " + (i + 1) + ": " + compression);
					if (!minimalSets.contains(compression)) {
//						System.err.println("Added Compression: " + compression);
						minimalSets.add(compression);
					}
				}
			}
			l++;
		}
		Set<NFE> result = new TreeSet<>();
		for (Set<NFE> minimalSet : minimalSets) {
			for (NFE t : minimalSet) {
				if (field.asVector(t).get(1).compareTo(Rationals.q().zero()) >= 0) {
					result.add(t);
				}
			}
		}
		return result;
	}

//	private Real logValue(Set<NFE> minimalSet, int index) {
//		Real result = null;
//		for (NFE alpha : minimalSet) {
//			Real value = field.logRepresentation(alpha).get(index);
//			if (result == null || value.compareTo(result) > 0) {
//				result = value;
//			}
//		}
//		return result;
//	}
//
//	private Vector<Real> logValues(Set<NFE> minimalSet) {
//		List<Real> result = new ArrayList<>();
//		for (int i = 0; i < field.realEmbeddings().size() + field.complexEmbeddings().size(); i++) {
//			result.add(logValue(minimalSet, i + 1));
//		}
//		return new Vector<>(result);
//	}

	private Real realValue(Set<NFE> minimalSet, EmbeddedNumberField<Real, Reals> embeddedField) {
		Real result = null;
		for (NFE alpha : minimalSet) {
			Real value = embeddedField.value(alpha);
			if (result == null || value.compareTo(result) > 0) {
				result = value;
			}
		}
		return result;
	}

	private Real complexValue(Set<NFE> minimalSet, EmbeddedNumberField<ComplexNumber, Complex> embeddedField) {
		Real result = null;
		Complex c = Complex.c(128);
		for (NFE alpha : minimalSet) {
			Real value = c.norm(embeddedField.embedding(alpha));
			if (result == null || value.compareTo(result) > 0) {
				result = value;
			}
		}
		return result;
	}

	private Real value(Set<NFE> minimalSet, int index) {
		if (index <= field.realEmbeddings().size()) {
			return realValue(minimalSet, field.realEmbeddings().get(index - 1));
		}
		return complexValue(minimalSet, field.complexEmbeddings().get(index - field.realEmbeddings().size() - 1));
	}

	private Vector<Real> values(Set<NFE> minimalSet) {
		List<Real> result = new ArrayList<>();
		for (int i = 0; i < field.realEmbeddings().size() + field.complexEmbeddings().size(); i++) {
			result.add(value(minimalSet, i + 1));
		}
		return new Vector<>(result);
	}

	private Real value(NFE alpha, int index) {
		if (index <= field.realEmbeddings().size()) {
			return field.realEmbeddings().get(index - 1).value(alpha);
		}
		EmbeddedNumberField<ComplexNumber, Complex> embedding = field.complexEmbeddings()
				.get(index - field.realEmbeddings().size() - 1);
		Complex c = Complex.c(128);
		return c.norm(embedding.embedding(alpha));
	}

	private Vector<Real> values(NFE alpha) {
		List<Real> result = new ArrayList<>();
		for (int i = 0; i < field.realEmbeddings().size() + field.complexEmbeddings().size(); i++) {
			result.add(value(alpha, i + 1));
		}
		return new Vector<>(result);
	}

	private Set<NFE> restriction(Set<NFE> minimalSet, int index) {
		Set<NFE> result = new TreeSet<>();
		Vector<Real> values = values(minimalSet);
		alphaLoop: for (NFE alpha : minimalSet) {
			Vector<Real> v = values(alpha);
			if (!v.get(index).equals(values.get(index))) {
				continue;
			}
			for (int i = 0; i < field.realEmbeddings().size() + field.complexEmbeddings().size(); i++) {
				if (i + 1 == index) {
					continue;
				}
				if (v.get(i + 1).compareTo(values.get(i + 1)) >= 0) {
					continue alphaLoop;
				}
			}
			result.add(alpha);
		}
		return result;
	}

	private List<Vector<Real>> cosineSineList(int number) {
		FiniteRealVectorSpace space = field.minkowskiEmbeddingSpace();
		Reals r = space.getValueField();
		List<Vector<Real>> result = new ArrayList<>();
		Real multiplier = r.divide(r.multiply(2, r.pi()), r.getInteger(number));
		for (int i = 0; i < number; i++) {
			Real angle = r.multiply(r.getInteger(i), multiplier);
			Vector<Real> cosineAndSine = r.cosineAndSine(angle);
			result.add(cosineAndSine);
		}
		return result;
	}

	private List<Dual<Real, Vector<Real>>> circle(int number, int quarterCircleIndex) {
		FiniteRealVectorSpace space = field.minkowskiEmbeddingSpace();
		Reals r = space.getValueField();
		List<Vector<Real>> cosineSineList = cosineSineList(number);
		int realEmbeddingIndex = field.realEmbeddings().size();
		List<Dual<Real, Vector<Real>>> constraints = new ArrayList<>();
		DualVectorSpace<Real, Vector<Real>> dual = space.getDual();
		for (int i = 0; i < space.dimension(); i++) {
			if (i < realEmbeddingIndex) {
				constraints.add(dual.getUnitVector(i + 1));
				constraints.add(dual.negative(dual.getUnitVector(i + 1)));
			} else {
				for (int j = 0; j < cosineSineList.size(); j++) {
					if (i + 1 == quarterCircleIndex && j >= number / 4 + 1) {
						break;
					}
					Vector<Real> cosineAndSine = cosineSineList.get(j);
					Real cosine = cosineAndSine.get(1);
					Real sine = cosineAndSine.get(2);
					List<Real> vector = new ArrayList<>();
					for (int k = 0; k < i; k++) {
						vector.add(r.zero());
					}
					vector.add(cosine);
					vector.add(sine);
					for (int k = i + 2; k < field.degree(); k++) {
						vector.add(r.zero());
					}
					constraints.add(space.getDual().canonicalIsomorphism(new Vector<>(vector)));
				}
				if (i + 1 == quarterCircleIndex) {
					constraints.add(dual.negative(dual.getUnitVector(i + 1)));
					constraints.add(dual.negative(dual.getUnitVector(i + 2)));
				}
				i++;
			}
		}
		return constraints;
	}

	private Set<NFE> expansion(Set<NFE> minimalSet, int index, FractionalIdeal ideal) {
		if (!restriction(minimalSet, index).isEmpty()) {
			return minimalSet;
		}
		Reals r = Reals.r(128);
		Real epsilon = r.power(r.getInteger(2), -16);
		Vector<Real> absoluteValues = values(minimalSet);
		int realEmbeddings = field.realEmbeddings().size();
		int complexEmbeddings = field.complexEmbeddings().size();
		List<Real> values = new ArrayList<>();
		int numComplexConstraints = 1 << (complexEmbeddings + 1);
		for (int i = 0; i < realEmbeddings + complexEmbeddings; i++) {
			if (i + 1 == index) {
				Real limit = r.abs(r.multiply(r.getInteger(100 * field.degree()), r.getInteger(field.discriminant())));
				if (i < realEmbeddings) {
					values.add(limit);
					values.add(r.negative(epsilon));
				} else {
					for (int j = 0; j < numComplexConstraints / 4 + 1; j++) {
						values.add(limit);
					}
					values.add(r.negative(epsilon));
					values.add(r.negative(epsilon));
				}
			} else {
				Real value = absoluteValues.get(i + 1);
				if (i >= realEmbeddings) {
					value = r.multiply(r.inverseSqrt(r.getInteger(2)), r.positiveSqrt(value));
				}
				value = r.subtract(value, epsilon);
				if (i < realEmbeddings) {
					values.add(value);
					values.add(value);
				} else {
					for (int j = 0; j < numComplexConstraints; j++) {
						values.add(value);
					}
				}
			}
		}
		int realIndex = index <= realEmbeddings ? index : realEmbeddings + 2 * (index - realEmbeddings) - 1;
		List<Dual<Real, Vector<Real>>> constraints = circle(numComplexConstraints, realIndex);
		DualVectorSpace<Real, Vector<Real>> dualVectorSpace = field.minkowskiEmbeddingSpace().getDual();
		Polytope<Real, Vector<Real>> polytope = new Polytope<>(field.minkowskiEmbeddingSpace(), constraints, values);
		Dual<Real, Vector<Real>> minimize = index <= realEmbeddings
				? dualVectorSpace.negative(dualVectorSpace.getUnitVector(realIndex))
				: dualVectorSpace.add(dualVectorSpace.negative(dualVectorSpace.getUnitVector(realIndex)),
						dualVectorSpace.negative(dualVectorSpace.getUnitVector(realIndex + 1)));
		NFE betterBorder = field.minkowskiEmbeddingSpace().latticeLinearProgram(polytope, minimize, ideal);
		constraints = circle(numComplexConstraints, -1);
		values.clear();
		for (int i = 0; i < realEmbeddings + complexEmbeddings; i++) {
			if (i + 1 == index) {
				Real limit = values(betterBorder).get(index);
				if (i >= realEmbeddings) {
					limit = r.positiveSqrt(limit);
					for (int j = 0; j < numComplexConstraints; j++) {
						values.add(limit);
					}
				} else {
					values.add(limit);
					values.add(limit);
				}
			} else {
				Real value = absoluteValues.get(i + 1);
				if (i >= realEmbeddings) {
					value = r.positiveSqrt(value);
					for (int j = 0; j < numComplexConstraints; j++) {
						values.add(value);
					}
				} else {
					values.add(value);
					values.add(value);
				}
			}
		}
		polytope = new Polytope<>(field.minkowskiEmbeddingSpace(), constraints, values);
		List<NFE> candidates = field.minkowskiEmbeddingSpace().latticePointsInPolytope(polytope, ideal);
		Real minValue = null;
		Set<NFE> add = new TreeSet<>();
		epsilon = r.getPowerOfTwo(-r.precision() / 2);
		candidateLoop: for (NFE candidate : candidates) {
			if (candidate.equals(zero())) {
				continue;
			}
			Vector<Real> candidateValues = values(candidate);
			for (int i = 0; i < realEmbeddings + complexEmbeddings; i++) {
				if (i + 1 == index) {
					continue;
				}
				if (candidateValues.get(i + 1).compareTo(r.subtract(absoluteValues.get(i + 1), epsilon)) >= 0) {
					continue candidateLoop;
				}
			}
			if (minValue == null || minValue.compareTo(candidateValues.get(index)) > 0) {
				add.clear();
				add.add(candidate);
				minValue = candidateValues.get(index);
			} else if (minValue.equals(candidateValues.get(index))) {
				add.add(candidate);
			}
		}
		Set<NFE> result = new TreeSet<>();
		result.addAll(minimalSet);
		result.addAll(add);
		return result;
	}

	private Set<NFE> compression(Set<NFE> minimalSet, int index) {
		Set<NFE> result = new TreeSet<>();
		Real value = value(minimalSet, index);
		for (NFE alpha : minimalSet) {
			Real v = value(alpha, index);
			if (v.compareTo(value) < 0) {
				result.add(alpha);
			}
		}
		return result;
	}

	@Override
	public FiniteField reduction(Ideal<NFE> maximalIdeal) {
		localize(maximalIdeal);
		return localizations.get(maximalIdeal).residueField();
	}

	@Override
	public FFE reduce(NFE t, Ideal<NFE> maximalIdeal) {
		return localize(maximalIdeal).reduce(t);
	}

	@Override
	public NFE lift(FFE s, Ideal<NFE> maximalIdeal) {
		return localize(maximalIdeal).lift(s);
	}

	public NFE centeredLift(FFE s, Ideal<NFE> maximalIdeal) {
		Integers z = Integers.z();
		IntE prime = ((NumberFieldIdeal) maximalIdeal).intGenerator();
		NFE lift = lift(s, maximalIdeal);
		Vector<IntE> asVector = asVector(lift);
		List<IntE> result = new ArrayList<>();
		for (int i = 0; i < field.degree(); i++) {
			result.add(z.centeredLift(z.reduce(asVector.get(i + 1), prime), prime));
		}
		return fromVector(new Vector<>(result));
	}

	@Override
	public LocalRingExtension<Fraction, PFE, NFE, LocalizedNumberField, PFE, FFE, FiniteField> localize(
			Ideal<NFE> maximalIdeal) {
		if (!localizations.containsKey(maximalIdeal)) {
			NumberFieldIdeal ideal = (NumberFieldIdeal) maximalIdeal;
			localizations.put(maximalIdeal,
					new LocalizedNumberField(field, ideal, Integers.z().localize(ideal.prime().getValue())));
		}
		return localizations.get(maximalIdeal).ringOfIntegers();
	}

	public LocalizedNumberField localizeAndQuotient(Ideal<NFE> maximalIdeal) {
		localize(maximalIdeal);
		return localizations.get(maximalIdeal);
	}

	@Override
	public NFE getEmbedding(IntE t) {
		return field.getEmbedding(t);
	}

	@Override
	public boolean isGeneratingAlgebra(List<NFE> s) {
		List<NFE> includingPowers = new ArrayList<>();
		includingPowers.add(one());
		for (NFE e : s) {
			NFE power = one();
			for (int i = 1; i < field.degree(); i++) {
				power = multiply(power, e);
				includingPowers.add(power);
			}
		}
		return isGeneratingModule(includingPowers);
	}

	@Override
	public List<NFE> getAlgebraGenerators() {
		return getModuleGenerators();
	}

	public class NumberFieldIdeal extends AbstractIdeal<NFE> implements Ideal<NFE> {
		private NFE uniformizer;
//		private NFE originalUniformizer;
		private IntE intGenerator;
		private Matrix<IntE> inBasis;
		private FreeSubModule<IntE, NFE> asSubModule;
		private MatrixModule<IntE> matrixModule;
		private List<NFE> generators;
		private OkutsuType<Fraction, PFE, PFE, FFE, FiniteField> type;
		private List<OkutsuType<Fraction, PFE, PFE, FFE, FiniteField>> allTypes;
		private int typeIndex;
		private List<List<UnivariatePolynomial<Fraction>>> integralBases;
		private boolean maximal;
		private FactorizationResult<Ideal<NFE>, Ideal<NFE>> idealFactorization;
		private NFE shortestElement;
		private boolean zeroIdeal;
		private boolean unitIdeal;
		private ModuloNumberFieldIdeal modOut;

		private NumberFieldIdeal(boolean zeroIdeal) {
			super(NumberFieldIntegers.this);
			if (zeroIdeal) {
				this.zeroIdeal = true;
				this.uniformizer = zero();
				this.intGenerator = Integers.z().zero();
				this.generators = Collections.emptyList();
				return;
			}
			this.unitIdeal = true;
			this.uniformizer = one();
			this.intGenerator = Integers.z().one();
			this.generators = Collections.singletonList(uniformizer);
			this.idealFactorization = new FactorizationResult<>(this, Collections.emptySortedMap());
		}

		private NumberFieldIdeal(
				TwoGeneratorIdeal<IntE, Fraction, PFE, PFE, FFE, FiniteField> ideal/*
																					 * DiscreteValuationRing<Fraction,
																					 * PFE> localized,
																					 * List<OkutsuType<Fraction, PFE,
																					 * PFE, FFE, FiniteField>> types,
																					 * int index,
																					 * List<List<UnivariatePolynomial<
																					 * Fraction>>> integral
																					 */) {
			super(NumberFieldIntegers.this);
			this.type = ideal.getType(); // types.get(index);
			this.allTypes = ideal.getTypes();// types;
			this.typeIndex = ideal.getIndex();// index;
			this.integralBases = ideal.getIntegralBases();// integral;
			this.uniformizer = field.fromPolynomial(ideal.getUniformizer());
			this.intGenerator = ideal.getIntGenerator();
			this.maximal = true;
//			Rationals q = Rationals.q();
//			UnivariatePolynomialRing<Fraction> polynomials = q.getUnivariatePolynomialRing();
//			UnivariatePolynomial<Fraction> uniformizer = null;
//			UnivariatePolynomial<Fraction> basisElement = null;
//			if (type.ramificationIndex() == 1) {
//				uniformizer = polynomials.getEmbedding(localized.uniformizer());
//				basisElement = integral.get(index).get(0);
//			} else {
//				boolean foundUniformizer = false;
//				boolean foundZero = false;
//				for (int i = 0; i < integral.get(index).size(); i++) {
//					Value value = type.valuation(integral.get(index).get(i));
//					if (!foundUniformizer && value.equals(Value.ONE)) {
//						foundUniformizer = true;
//						uniformizer = integral.get(index).get(i);
//					} else if (!foundZero && value.equals(Value.ZERO)) {
//						foundZero = true;
//						basisElement = integral.get(index).get(i);
//					}
//					if (foundZero && foundUniformizer) {
//						break;
//					}
//				}
//				if (!foundUniformizer) {
//					System.err.println("Did not find uniformizer, but kinda expected to find one!");
//					uniformizer = type.lift(type.reduction().extension().one(), 1);
//				}
//				if (!foundZero) {
//					throw new ArithmeticException("Did not find zero valuation element in integral basis!");
//				}
//			}
//			UnivariatePolynomial<Fraction> generator = polynomials.multiply(uniformizer, basisElement);
//			for (int i = 0; i < types.size(); i++) {
//				if (i == index) {
//					continue;
//				}
//				OkutsuType<Fraction, PFE, PFE, FFE, FiniteField> otherType = types.get(i);
//				if (!otherType.valuation(generator).equals(Value.ZERO)) {
//					for (int j = 0; j < integral.get(i).size(); j++) {
//						if (otherType.valuation(integral.get(i).get(j)).equals(Value.ZERO)) {
//							generator = polynomials.add(generator, integral.get(i).get(j));
//							break;
//						}
//					}
//				}
//			}
//			for (int i = 0; i < types.size(); i++) {
//				if (i == index) {
//					if (!types.get(i).valuation(generator).equals(Value.ONE)) {
//						throw new ArithmeticException("Not a uniformizer of the prime ideal!");
//					}
//				} else {
//					if (!types.get(i).valuation(generator).equals(Value.ZERO)) {
//						throw new ArithmeticException("Not a unit mod the other prime ideals!");
//					}
//				}
//			}
//			this.uniformizer = field.fromPolynomial(generator);
//			this.intGenerator = localized.uniformizer().asInteger();
			this.idealFactorization = new FactorizationResult<>(getUnitIdeal(), SingletonSortedMap.map(this, 1));
			init();
			for (int i = 0; i < allTypes.size(); i++) {
				if (i == typeIndex) {
					if (!allTypes.get(i).valuation(this.uniformizer.asPolynomial()).equals(Value.ONE)) {
						throw new ArithmeticException("Not a uniformizer of the prime ideal!");
					}
				} else {
					if (!allTypes.get(i).valuation(this.uniformizer.asPolynomial()).equals(Value.ZERO)) {
						throw new ArithmeticException("Not a unit mod the other prime ideals!");
					}
				}
			}
		}

		private NumberFieldIdeal(Map<IntE, FactorizationResult<Ideal<NFE>, Ideal<NFE>>> factorization) {
			super(NumberFieldIntegers.this);
			SortedMap<Ideal<NFE>, Integer> factors = new TreeMap<>();
			Map<IntE, Integer> hMap = new TreeMap<>();
			Map<IntE, NFE> alphaMap = new TreeMap<>();
			for (IntE prime : factorization.keySet()) {
				FactorizationResult<Ideal<NFE>, Ideal<NFE>> primeFactorization = factorization.get(prime);
				NFE alpha = one();
				int h = 0;
				for (Ideal<NFE> primeIdeal : primeFactorization.primeFactors()) {
					NumberFieldIdeal ideal = (NumberFieldIdeal) primeIdeal;
					if (!prime.equals(ideal.intGenerator)) {
						throw new ArithmeticException("Prime does match prime ideal!");
					}
					if (!ideal.maximal) {
						throw new ArithmeticException("Prime Ideal factor not prime!");
					}
					int multiplicity = primeFactorization.multiplicity(primeIdeal);
					factors.put(primeIdeal, multiplicity);
					alpha = multiply(alpha, power(ideal.uniformizer, multiplicity));
					h = Math.max(h, MiscAlgorithms.DivRoundUp(multiplicity, ideal.type.ramificationIndex()));
				}
				hMap.put(prime, h);
				alphaMap.put(prime, alpha);
			}
			Integers z = Integers.z();
			this.intGenerator = z.one();
			this.uniformizer = zero();
			for (IntE prime : factorization.keySet()) {
				this.intGenerator = z.multiply(intGenerator, z.power(prime, hMap.get(prime)));
				NFE alpha = alphaMap.get(prime);
				for (IntE otherPrime : factorization.keySet()) {
					if (prime.equals(otherPrime)) {
						continue;
					}
					alpha = multiply(z.power(otherPrime, hMap.get(otherPrime) + 1), alpha);
				}
				this.uniformizer = add(this.uniformizer, alpha);
			}
			this.idealFactorization = new FactorizationResult<>(getUnitIdeal(), factors);
			this.maximal = false;
			if (idealFactorization.isIrreducible()) {
				throw new ArithmeticException("Prime Ideals need to use different constructor!");
			}
			init();
		}

		private boolean proposeNewUniformizer(NFE proposedUniformizer) {
			if (!asSubModule.contains(proposedUniformizer)) {
				return false;
			}
			if (proposedUniformizer.asPolynomial().degree() <= 0) {
				return false;
			}
			if (field.norm(proposedUniformizer).asInteger().getValue().abs()
					.compareTo(field.norm(uniformizer).asInteger().getValue().abs()) >= 0) {
				return false;
			}
			Integers z = Integers.z();
			if (maximal) {
				if (!type.valuation(proposedUniformizer.asPolynomial()).equals(Value.ONE)) {
					if (type.ramificationIndex() != 1) {
						return false;
					}
					boolean success1 = proposeNewUniformizer(add(proposedUniformizer, getInteger(intGenerator)));
					boolean success2 = proposeNewUniformizer(subtract(proposedUniformizer, getInteger(intGenerator)));
					return success1 || success2;
				}
				for (int i = 0; i < allTypes.size(); i++) {
					if (i == typeIndex) {
						continue;
					}
					OkutsuType<Fraction, PFE, PFE, FFE, FiniteField> otherType = allTypes.get(i);
					if (!otherType.valuation(proposedUniformizer.asPolynomial()).equals(Value.ZERO)) {
						for (int j = 0; j < integralBases.get(i).size(); j++) {
							if (otherType.valuation(integralBases.get(i).get(j)).equals(Value.ZERO)) {
								proposedUniformizer = field.add(proposedUniformizer,
										field.fromSmallDegreePolynomial(integralBases.get(i).get(j)));
								break;
							}
						}
						return proposeNewUniformizer(proposedUniformizer);
					}
				}
			} else {
				for (IntE prime : z.uniqueFactorization(intGenerator).primeFactors()) {
					for (NumberFieldIdeal primeIdeal : idealsOver(prime)) {
						if (!primeIdeal.maximumPowerContains(proposedUniformizer)
								.equals(new Value(idealFactorization.multiplicity(primeIdeal)))) {
							return false;
						}
					}
				}
			}
			List<NFE> asList = new ArrayList<>();
			for (NFE b : NumberFieldIntegers.this.getModuleGenerators()) {
				NFE multiplied = multiply(b, proposedUniformizer);
				asList.add(multiplied);
			}
			for (NFE b : NumberFieldIntegers.this.getModuleGenerators()) {
				NFE multiplied = multiply(intGenerator, b);
				asList.add(multiplied);
			}
			FreeSubModule<IntE, NFE> proposedAsSubModule = new FreeSubModule<>(NumberFieldIntegers.this, asList);
			if (!proposedAsSubModule.contains(uniformizer)) {
				return false;
			}
			this.uniformizer = proposedUniformizer;
			this.asSubModule = proposedAsSubModule;
			return true;
		}

		private void init() {
			if (this.uniformizer.equals(zero())) {
				this.uniformizer = field.getEmbedding(intGenerator);
			}
			List<NFE> asList = new ArrayList<>();
			for (NFE b : NumberFieldIntegers.this.getModuleGenerators()) {
				NFE multiplied = multiply(b, uniformizer);
				asList.add(multiplied);
			}
			for (NFE b : NumberFieldIntegers.this.getModuleGenerators()) {
				NFE multiplied = multiply(intGenerator, b);
				asList.add(multiplied);
			}
			matrixModule = new MatrixModule<>(Integers.z(), field.degree(), 2 * field.degree());
			asSubModule = new FreeSubModule<>(NumberFieldIntegers.this, asList);
			if (!asSubModule.contains(uniformizer)) {
				throw new ArithmeticException("Uniformizer not in submodule?");
			}
			if (!asSubModule.contains(field.getEmbedding(intGenerator))) {
				throw new ArithmeticException("int generator not in submodule?");
			}
			boolean skipLLL = field.degree() > 128;
			if (isElement(field.divide(uniformizer, field.getEmbedding(intGenerator)))
					|| isElement(field.divide(field.getEmbedding(intGenerator), uniformizer))) {
				skipLLL = true;
			}
//			for (NFE b : asSubModule.getBasis()) {
//				skipLLL = skipLLL || field.norm(b).getNumerator().getValue().bitLength() > 65536;
//			}
			List<NFE> basis = asSubModule.getBasis();
			if (!skipLLL) {
				basis = sublatticeReduction(basis);
			}
			if (isElement(field.divide(uniformizer, field.getEmbedding(intGenerator)))) {
				this.uniformizer = field.getEmbedding(intGenerator);
			} else {
				for (NFE b : basis) {
					proposeNewUniformizer(b);
				}
			}
			if (isElement(field.divide(uniformizer, field.getEmbedding(intGenerator)))) {
				this.uniformizer = field.getEmbedding(intGenerator);
			}
			List<Vector<IntE>> asVectorList = new ArrayList<>();
			for (NFE b : NumberFieldIntegers.this.getModuleGenerators()) {
				NFE multiplied = multiply(b, uniformizer);
				asVectorList.add(NumberFieldIntegers.this.asVector(multiplied));
			}
			for (NFE b : NumberFieldIntegers.this.getModuleGenerators()) {
				NFE multiplied = multiply(intGenerator, b);
				asVectorList.add(NumberFieldIntegers.this.asVector(multiplied));
			}
			this.inBasis = Matrix.fromColumns(asVectorList);
			this.generators = new ArrayList<>();
			this.generators.add(uniformizer);
			if (!isElement(field.divide(field.getEmbedding(intGenerator), uniformizer))) {
				this.generators.add(field.getEmbedding(intGenerator));
			}
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
		public boolean isPrincipal() {
			if (generators.size() == 1) {
				return true;
			}
			NumberFieldIdeal representative = reducedRepresentative(this).clearDenominator().getFirst();
			for (NumberFieldIdeal ideal : field.idealClassGroup().neutral().alternativeRepresentatives()) {
				if (ideal.equals(representative)) {
					return true;
				}
			}
			return false;
		}

		public ModuloNumberFieldIdeal modOut() {
			if (modOut == null) {
				this.modOut = new ModuloNumberFieldIdeal(NumberFieldIntegers.this, this);
			}
			return modOut;
		}

		public NFE shortestElement() {
			if (shortestElement == null) {
				shortestElement = field.minkowskiEmbeddingSpace()
						.latticeReduction(field.idealGroup().getEmbedding(this)).get(0);
			}
			return shortestElement;
		}

		@Override
		public NFE principalGenerator() {
			if (generators.size() == 1) {
				return generators.get(0);
			}
			FractionalIdeal reducedRepresentative = reducedRepresentative(this);
			NumberFieldIdeal representative = reducedRepresentative.clearDenominator().getFirst();
			List<NumberFieldIdeal> alternatives = field.idealClassGroup().neutral().alternativeRepresentatives();
			for (int i = 0; i < alternatives.size(); i++) {
				NumberFieldIdeal ideal = alternatives.get(i);
				if (ideal.equals(representative)) {
					return field.multiply(unitIdealMultipliers.get(i), shortestElement());
				}
			}
			throw new ArithmeticException("Not a principal ideal!");
		}

		public FreeSubModule<IntE, NFE> asSubModule() {
			if (unitIdeal) {
				return new FreeSubModule<>(NumberFieldIntegers.this, NumberFieldIntegers.this.getModuleGenerators());
			} else if (zeroIdeal) {
				return new FreeSubModule<>(NumberFieldIntegers.this, Collections.emptyList());
			}
			return asSubModule;
		}

		@Override
		public boolean isPrimary() {
			return isPrime() || idealFactorization.primeFactors().size() == 1;
		}

		@Override
		public boolean isPrime() {
			return zeroIdeal || isMaximal();
		}

		@Override
		public boolean isMaximal() {
			return maximal;
		}

		public IntE intGenerator() {
			return intGenerator;
		}

		public NFE uniformizer() {
			return uniformizer;
		}

//		@Override
//		public List<Vector<NFE>> nonTrivialCombinations(List<NFE> s) {
//			List<Vector<IntE>> asVectorList = new ArrayList<>();
//			for (NFE e : s) {
//				for (NFE basisVector : NumberFieldIntegers.this.getModuleGenerators()) {
//					asVectorList.add(NumberFieldIntegers.this.asVector(multiply(e, basisVector)));
//				}
//			}
//			List<List<IntE>> nonTrivialCombinations = asFreeModule.nonTrivialCombinations(asVectorList);
//			List<List<NFE>> result = new ArrayList<>();
//			int degree = field.degree();
//			for (List<IntE> combination : nonTrivialCombinations) {
//				List<NFE> row = new ArrayList<>();
//				for (int i = 0; i < s.size(); i++) {
//					row.add(NumberFieldIntegers.this
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
		public List<NFE> generators() {
			return generators;
		}

		@Override
		public List<NFE> generate(NFE t) {
			if (zeroIdeal) {
				return Collections.emptyList();
			}
			NFE inIdeal = field.subtract(t, residue(t));
			if (generators.size() == 1) {
				return Collections.singletonList(field.divide(inIdeal, uniformizer));
			}
			Vector<IntE> inIntegerBasis = NumberFieldIntegers.this.asVector(inIdeal);
			Vector<IntE> inIdealBasis = matrixModule.solve(inBasis, inIntegerBasis);
			List<NFE> g = new ArrayList<>();
			g.add(NumberFieldIntegers.this.fromVector(new Vector<>(inIdealBasis.asList().subList(0, field.degree()))));
			g.add(NumberFieldIntegers.this
					.fromVector(new Vector<>(inIdealBasis.asList().subList(field.degree(), 2 * field.degree()))));
			return g;
		}

		public NFE round(NFE t, int accuracy) {
			Value neededValue = new Value(accuracy);
			NFE rounded = zero();
			while (true) {
				Value value = type.valuation(t.asPolynomial());
				if (value.compareTo(neededValue) >= 0) {
					return rounded;
				}
				FFE reduced = type.reduce(t.asPolynomial());
				NFE lifted = field.fromPolynomial(type.lift(reduced, value.value()));
				rounded = add(rounded, lifted);
				t = subtract(t, lifted);
			}
		}

		@Override
		public NFE residue(NFE t) {
			if (zeroIdeal) {
				return t;
			}
			if (unitIdeal) {
				return zero();
			}
			return modOut().lift(modOut().reduce(t));
		}

		@Override
		public boolean contains(NFE t) {
			if (t.equals(zero()) || unitIdeal) {
				return true;
			}
			if (zeroIdeal) {
				return false;
			}
			for (Ideal<NFE> primeIdeal : idealFactorization.primeFactors()) {
				NumberFieldIdeal ideal = (NumberFieldIdeal) primeIdeal;
				Value value = ideal.type.valuation(t.asPolynomial());
				if (value.compareTo(new Value(idealFactorization.multiplicity(primeIdeal))) < 0) {
					return false;
				}
			}
			return true;
		}

		@Override
		public Value maximumPowerContains(NFE t) {
			if (t.equals(zero()) || unitIdeal) {
				return Value.INFINITY;
			}
			if (zeroIdeal) {
				return Value.ZERO;
			}
			Value value = Value.INFINITY;
			for (Ideal<NFE> primeIdeal : idealFactorization.primeFactors()) {
				NumberFieldIdeal ideal = (NumberFieldIdeal) primeIdeal;
				value = value.min(new Value(
						ideal.type.valuation(t.asPolynomial()).value() / idealFactorization.multiplicity(primeIdeal)));
			}
			return value;
		}

		@Override
		public boolean contains(Ideal<NFE> other) {
			NumberFieldIdeal otherIdeal = (NumberFieldIdeal) other;
			if (otherIdeal.zeroIdeal) {
				return true;
			}
			if (zeroIdeal) {
				return false;
			}
			if (unitIdeal) {
				return true;
			}
			if (otherIdeal.unitIdeal) {
				return false;
			}
			if (otherIdeal.maximal && !maximal) {
				return false;
			}
			if (otherIdeal.maximal && otherIdeal.uniformizer.equals(uniformizer)
					&& otherIdeal.intGenerator.equals(intGenerator)) {
				return true;
			}
			if (otherIdeal.maximal) {
				return false;
			}
			for (Ideal<NFE> otherPrimeFactor : otherIdeal.idealFactorization.primeFactors()) {
				if (idealFactorization.multiplicity(otherPrimeFactor) < otherIdeal.idealFactorization
						.multiplicity(otherPrimeFactor)) {
					return false;
				}
			}
			return true;
		}

		@Override
		public int compareTo(Ideal<NFE> other) {
			NumberFieldIdeal otherIdeal = (NumberFieldIdeal) other;
			if (!maximal || !otherIdeal.maximal) {
				Set<Ideal<NFE>> ideals = new TreeSet<>();
				FactorizationResult<Ideal<NFE>, Ideal<NFE>> thisFactorization = idealFactorization(this);
				FactorizationResult<Ideal<NFE>, Ideal<NFE>> otherFactorization = idealFactorization(other);
				ideals.addAll(thisFactorization.primeFactors());
				ideals.addAll(otherFactorization.primeFactors());
				for (Ideal<NFE> ideal : ideals) {
					int cmp = thisFactorization.multiplicity(ideal) - otherFactorization.multiplicity(ideal);
					if (cmp != 0) {
						return cmp;
					}
				}
				return 0;
			}
			int cmp = intGenerator.compareTo(otherIdeal.intGenerator);
			if (cmp != 0) {
				return cmp;
			}
			return uniformizer.compareTo(otherIdeal.uniformizer);
		}

		@Override
		public boolean equalsIdeal(Ideal<NFE> other) {
			NumberFieldIdeal otherIdeal = (NumberFieldIdeal) other;
			if (!maximal || !otherIdeal.maximal) {
				return super.equalsIdeal(other);
			}
			return compareTo(other) == 0;
		}

		public IntE prime() {
			if (!maximal) {
				throw new ArithmeticException("not a maximal ideal!");
			}
			return intGenerator;
		}

		public OkutsuType<Fraction, PFE, PFE, FFE, FiniteField> type() {
			if (!maximal) {
				throw new ArithmeticException("not a maximal ideal!");
			}
			return type;
		}

		public IntE norm() {
			Integers z = Integers.z();
			if (maximal) {
				return z.getInteger(type.reduction().extension().getNumberOfElements());
			}
			IntE norm = z.one();
			for (Ideal<NFE> primeIdeal : idealFactorization.primeFactors()) {
				NumberFieldIdeal ideal = (NumberFieldIdeal) primeIdeal;
				norm = z.multiply(norm, z.power(ideal.norm(), idealFactorization.multiplicity(primeIdeal)));
			}
			return norm;
		}

		public NumberFieldIdeal normOver(FieldEmbedding<Fraction, NFE, NumberField> fieldEmbedding) {
			SortedMap<NumberFieldIdeal, Integer> factorization = new TreeMap<>();
			for (Ideal<NFE> primeFactor : idealFactorization(this).primeFactors()) {
				NumberFieldIdeal lowerPrimeIdeal = intersectToLowerField(primeFactor, fieldEmbedding);
				factorization.put(lowerPrimeIdeal, ((NumberFieldIdeal) primeFactor).type().residueDegree()
						/ lowerPrimeIdeal.type().residueDegree() * idealFactorization(this).multiplicity(primeFactor));
			}
			return fieldEmbedding.getEmbeddedField().maximalOrder().fromFactorization(factorization);
		}

		public IntE discriminant() {
			return NumberFieldIntegers.this.discriminant(asSubModule.getBasis());
		}

		public Real minkowskiBound() {
			Reals r = Reals.r(1024);
			Real twoByPiToS = r.power(r.divide(r.getInteger(2), r.pi()), field.complexEmbeddings().size());
			Real sqrtDiscriminant = r.positiveSqrt(r.abs(r.getEmbedding(field.discriminant())));
			Real norm = r.getEmbedding(norm());
			return r.multiply(twoByPiToS, sqrtDiscriminant, norm);
		}
	}
}
