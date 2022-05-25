package fields.numberfields;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import fields.exceptions.InfinityException;
import fields.helper.AbstractElement;
import fields.helper.AbstractModule;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.interfaces.Group;
import fields.interfaces.Ideal;
import fields.interfaces.Module;
import fields.numberfields.PicardGroup.IdealClass;
import fields.numberfields.NumberField.NFE;
import fields.numberfields.NumberFieldIntegers.NumberFieldIdeal;
import fields.vectors.FreeModule;
import fields.vectors.GenericPIDModule;
import fields.vectors.GenericPIDModule.Mod;
import fields.vectors.Vector;

public class PicardGroup extends AbstractModule<IntE, IdealClass>
		implements Group<IdealClass>, Module<IntE, IdealClass> {
	private NumberFieldIntegers order;
	private IdealGroup idealGroup;
	private Map<NumberFieldIdeal, IdealClass> reducedMapping;
	private int nextIndex;
	private List<IdealClass> generators;
	private List<Integer> orders;
	private IdealClass neutral;
	private Map<Integer, Map<Integer, IdealClass>> operationMap;
	private Map<Integer, IdealClass> inverseMap;
	private GenericPIDModule<IntE, Vector<IntE>> asGenericModule;

	public static class IdealClass extends AbstractElement<IdealClass> {
		private NumberFieldIdeal representative;
		private List<NumberFieldIdeal> alternativeRepresentatives;
		private int index;
		private Vector<IntE> asVector;

		private IdealClass(int index, NumberFieldIdeal representative,
				List<NumberFieldIdeal> alternativeRepresentatives) {
			this.index = index;
			this.representative = representative;
			this.alternativeRepresentatives = alternativeRepresentatives;
		}

		@Override
		public int compareTo(IdealClass o) {
			return index - o.index;
		}

		@Override
		public String toString() {
			return "[" + representative + "]";
		}

		public NumberFieldIdeal representative() {
			return representative;
		}

		List<NumberFieldIdeal> alternativeRepresentatives() {
			return alternativeRepresentatives;
		}

		public boolean isPrincipal() {
			return representative.isPrincipal();
		}
	}

	PicardGroup(NumberFieldIntegers order) {
		this.order = order;
		this.idealGroup = order.numberField().idealGroup();
		this.reducedMapping = new TreeMap<>();
		this.nextIndex = 0;
		this.neutral = getEmbedding(idealGroup.neutral());
		this.operationMap = new TreeMap<>();
		this.inverseMap = new TreeMap<>();
	}

	public IdealClass getEmbedding(FractionalIdeal ideal) {
		return getEmbedding(ideal.clearDenominator().getFirst());
	}

	public IdealClass getEmbedding(Ideal<NFE> ideal) {
		NumberFieldIdeal representative = order.reducedRepresentative(ideal).clearDenominator().getFirst();
		if (!reducedMapping.containsKey(representative)) {
			List<NumberFieldIdeal> representatives = order.idealEquivalenceClass(ideal);
			IdealClass idealClass = new IdealClass(nextIndex, representative, representatives);
			nextIndex++;
			for (NumberFieldIdeal alternate : representatives) {
				reducedMapping.put(alternate, idealClass);
			}
		}
		return reducedMapping.get(representative);

	}

	private void computeAllElements() {
		if (generators != null) {
			return;
		}
		Integers z = Integers.z();
		List<List<IntE>> syzygies = new ArrayList<>();
		List<IdealClass> unoptimizedGenerators = unoptimizedGenerators();
		FreeModule<IntE> free = new FreeModule<IntE>(z, unoptimizedGenerators.size());
		neutral().asVector = free.zero();
		List<IdealClass> generatedElements = new ArrayList<>();
		generatedElements.add(neutral());
		for (int i = 0; i < unoptimizedGenerators.size(); i++) {
			IdealClass generator = unoptimizedGenerators.get(i);
			IdealClass generatorPower = generator;
			int power = 1;
			while (!generatorPower.equals(neutral())) {
				Vector<IntE> powerVector = free.scalarMultiply(power, free.getUnitVector(i + 1));
				if (generatorPower.asVector == null) {
					List<IdealClass> elementsCopy = new ArrayList<>();
					elementsCopy.addAll(generatedElements);
					for (IdealClass element : elementsCopy) {
						IdealClass cl = operate(element, generatorPower);
						if (cl.asVector != null) {
							continue;
						}
						generatedElements.add(cl);
						cl.asVector = free.add(element.asVector, powerVector);
					}
				} else {
					Vector<IntE> syzygy = free.subtract(generatorPower.asVector, powerVector);
					syzygies.add(syzygy.asList());
				}
				generatorPower = operate(generatorPower, generator);
				power++;
			}
			syzygies.add(free.scalarMultiply(power, free.getUnitVector(i + 1)).asList());
		}
		this.asGenericModule = GenericPIDModule.fromSyzygies(free, syzygies);
		this.generators = new ArrayList<>();
		this.orders = new ArrayList<>();
		syzygies = new ArrayList<>();
		free = null;
		for (int i = 0; i < asGenericModule.diagonalBasis().size(); i++) {
			if (z.isUnit(asGenericModule.diagonalRanks().get(i))) {
				continue;
			}
			if (free == null) {
				free = new FreeModule<>(z, asGenericModule.diagonalBasis().size() - i);
			}
			Vector<IntE> basisVector = asGenericModule.lift(asGenericModule.diagonalBasis().get(i));
			IdealClass generator = neutral();
			for (int j = 0; j < unoptimizedGenerators.size(); j++) {
				generator = operate(power(basisVector.get(j+1).intValueExact(), unoptimizedGenerators.get(j)), generator);
			}
			generators.add(generator);
			orders.add(Math.abs(asGenericModule.diagonalRanks().get(i).intValueExact()));
			syzygies.add(free
					.scalarMultiply(orders.get(free.dimension() - asGenericModule.diagonalBasis().size() + i),
							free.getUnitVector(free.dimension() - asGenericModule.diagonalBasis().size() + i + 1))
					.asList());
		}
		if (free == null) {
			return;
		}
		asGenericModule = GenericPIDModule.fromSyzygies(free, syzygies);
	}

	@Override
	public String toString() {
		return "Cl(" + order.numberField() + ")";
	}

	@Override
	public Integers getRing() {
		return Integers.z();
	}

	@Override
	public Exactness exactness() {
		return order.exactness();
	}

	public GenericPIDModule<IntE, Vector<IntE>> asIntModule() {
		computeAllElements();
		return asGenericModule;
	}

	@Override
	public IdealClass getRandomElement() {
		return fromVector(asIntModule().lift(asIntModule().getRandomElement()));
	}

	@Override
	public boolean isFinite() {
		return true;
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		return asIntModule().getNumberOfElements();
	}

	@Override
	public Iterator<IdealClass> iterator() {
		return new Iterator<>() {
			private Iterator<Mod<Vector<IntE>>> it = asIntModule().iterator();

			@Override
			public boolean hasNext() {
				return it.hasNext();
			}

			@Override
			public IdealClass next() {
				return fromVector(asIntModule().lift(it.next()));
			}
		};
	}

	@Override
	public IdealClass neutral() {
		return neutral;
	}

	@Override
	public IdealClass zero() {
		return neutral();
	}

	@Override
	public IdealClass inverse(IdealClass t) {
		if (!inverseMap.containsKey(t.index)) {
			IdealClass cl = getEmbedding(idealGroup.inverse(idealGroup.getEmbedding(t.representative)));
			inverseMap.put(t.index, cl);
		}
		return inverseMap.get(t.index);
	}

	@Override
	public IdealClass negative(IdealClass s) {
		return inverse(s);
	}

	@Override
	public IdealClass operate(IdealClass t1, IdealClass t2) {
		if (!this.operationMap.containsKey(t1.index)) {
			this.operationMap.put(t1.index, new TreeMap<>());
		}
		if (!this.operationMap.containsKey(t2.index)) {
			this.operationMap.put(t2.index, new TreeMap<>());
		}
		if (!this.operationMap.get(t1.index).containsKey(t2.index)) {
			IdealClass cl = getEmbedding(order.multiply(t1.representative, t2.representative));
			this.operationMap.get(t1.index).put(t2.index, cl);
			this.operationMap.get(t2.index).put(t1.index, cl);
		}
		return operationMap.get(t1.index).get(t2.index);
	}

	@Override
	public IdealClass add(IdealClass s1, IdealClass s2) {
		return operate(s1, s2);
	}

	@Override
	public IdealClass scalarMultiply(IntE t, IdealClass s) {
		return power(t.getValue(), s);
	}

	@Override
	public boolean isFree() {
		return false;
	}

	@Override
	public Ideal<IntE> annihilator() {
		return asGenericModule.annihilator();
	}

	@Override
	public boolean isLinearIndependent(List<IdealClass> s) {
		return s.isEmpty();
	}

	private Mod<Vector<IntE>> asIntModuleElement(IdealClass s) {
		return asIntModule().reduce(asVector(s));
	}

	private List<Mod<Vector<IntE>>> asIntModuleElementList(List<IdealClass> s) {
		List<Mod<Vector<IntE>>> result = new ArrayList<>();
		for (IdealClass c : s) {
			result.add(asIntModuleElement(c));
		}
		return result;
	}

	@Override
	public boolean isGeneratingModule(List<IdealClass> s) {
		return asGenericModule.isGeneratingModule(asIntModuleElementList(s));
	}

	@Override
	public List<List<IntE>> nonTrivialCombinations(List<IdealClass> s) {
		return asIntModule().nonTrivialCombinations(asIntModuleElementList(s));
	}

	private List<IdealClass> unoptimizedGenerators() {
		Integers z = Integers.z();
		List<IdealClass> generators = new ArrayList<>();
		IntE minkowskiLimit = order.numberField().minkowskiBound().roundDown();
		for (IntE prime : z.setOfPrimes()) {
			if (prime.compareTo(minkowskiLimit) > 0) {
				break;
			}
			for (NumberFieldIdeal ideal : order.idealsOver(prime)) {
				IdealClass generator = getEmbedding(ideal);
				if (generator.equals(neutral())) {
					continue;
				}
				generators.add(generator);
			}
		}
		if (generators.isEmpty()) {
			generators.add(neutral());
		}
		return generators;
	}

	@Override
	public List<IdealClass> getModuleGenerators() {
		if (generators == null) {
			computeAllElements();
		}
		return generators;
	}

	public List<Integer> getOrders() {
		 asIntModule();
		 return orders;
	}

	@Override
	public List<List<IntE>> getModuleGeneratorRelations() {
		return asIntModule().getModuleGeneratorRelations();
	}

	@Override
	public Vector<IntE> asVector(IdealClass s) {
		computeAllElements();
		return s.asVector;
	}
}
