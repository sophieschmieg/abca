package fields.numberfields;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import fields.exceptions.InfinityException;
import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.helper.AbstractElement;
import fields.helper.AbstractModule;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.interfaces.Group;
import fields.interfaces.Ideal;
import fields.interfaces.Module;
import fields.numberfields.NumberField.NFE;
import fields.numberfields.NumberFieldIntegers.NumberFieldIdeal;
import fields.numberfields.NumberFieldOrder.NumberFieldOrderIdeal;
import fields.numberfields.PicardGroup.OrderIdealClass;
import fields.vectors.FreeModule;
import fields.vectors.GenericPIDModule;
import fields.vectors.GenericPIDModule.Mod;
import fields.vectors.Vector;

public class PicardGroup extends AbstractModule<IntE, OrderIdealClass>
		implements Group<OrderIdealClass>, Module<IntE, OrderIdealClass> {
	private NumberFieldOrder order;
	private OrderIdealGroup idealGroup;
	private Map<FractionalOrderIdeal, OrderIdealClass> reducedMapping;
	private List<OrderIdealClass> elements;
	private List<OrderIdealClass> generators;
	private List<Integer> orders;
	private OrderIdealClass neutral;
	private Map<Integer, Map<Integer, OrderIdealClass>> operationMap;
	private Map<Integer, OrderIdealClass> inverseMap;
	private GenericPIDModule<IntE, Vector<IntE>> asGenericModule;

	public static class OrderIdealClass extends AbstractElement<OrderIdealClass> {
		private FractionalOrderIdeal representative;
		private int index;
		private Vector<IntE> asVector;

		private OrderIdealClass(int index, FractionalOrderIdeal representative) {
			this.index = index;
			this.representative = representative;
		}

		@Override
		public int compareTo(OrderIdealClass o) {
			return index - o.index;
		}

		@Override
		public String toString() {
			return "[" + representative + "]";
		}

		public FractionalOrderIdeal representative() {
			return representative;
		}

		public boolean isPrincipal() {
			return representative.isPrincipal();
		}
	}

	PicardGroup(NumberFieldOrder order) {
		this.order = order;
		this.idealGroup = order.idealGroup();
		this.elements = new ArrayList<>();
		this.reducedMapping = new TreeMap<>();
		this.neutral = getEmbedding(idealGroup.neutral());
		this.operationMap = new TreeMap<>();
		this.inverseMap = new TreeMap<>();
	}

	public OrderIdealClass getEmbedding(Ideal<NFE> ideal) {
		return getEmbedding(idealGroup.getEmbedding(ideal));
	}

	public OrderIdealClass getEmbedding(FractionalOrderIdeal ideal) {
		NFE minimal = ideal.getVectorSpace().latticeReduction(ideal).get(0);
		FractionalOrderIdeal representative = idealGroup.operate(ideal,
				idealGroup.getPrincipalIdeal(order.numberField().inverse(minimal)));
		if (!reducedMapping.containsKey(representative)) {
			FractionalOrderIdeal inverse = idealGroup.inverse(representative);
			for (OrderIdealClass element : elements) {
				FractionalOrderIdeal divided = idealGroup.operate(element.representative(), inverse);
				if (divided.isPrincipal()) {
					reducedMapping.put(representative, element);
					return element;
				}
			}
			OrderIdealClass idealClass = new OrderIdealClass(elements.size(), representative);
			elements.add(idealClass);
			reducedMapping.put(representative, idealClass);
		}
		return reducedMapping.get(representative);

	}

	private void computeAllElements() {
		if (generators != null) {
			return;
		}
		Integers z = Integers.z();
		List<List<IntE>> syzygies = new ArrayList<>();
		List<OrderIdealClass> unoptimizedGenerators = unoptimizedGenerators();
		FreeModule<IntE> free = new FreeModule<IntE>(z, unoptimizedGenerators.size());
		neutral().asVector = free.zero();
		for (int i = 0; i < unoptimizedGenerators.size(); i++) {
			OrderIdealClass generator = unoptimizedGenerators.get(i);
			OrderIdealClass generatorPower = generator;
			int power = 1;
			while (!generatorPower.equals(neutral())) {
				Vector<IntE> powerVector = free.scalarMultiply(power, free.getUnitVector(i + 1));
				if (generatorPower.asVector == null) {
					List<OrderIdealClass> elementsCopy = new ArrayList<>();
					elementsCopy.addAll(elements);
					for (OrderIdealClass element : elementsCopy) {
						OrderIdealClass cl = operate(element, generatorPower);
						if (cl.asVector != null) {
							continue;
						}
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
			OrderIdealClass generator = neutral();
			for (int j = 0; j < unoptimizedGenerators.size(); j++) {
				generator = operate(power(basisVector.get(j + 1).intValueExact(), unoptimizedGenerators.get(j)),
						generator);
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
		GenericPIDModule<IntE, Vector<IntE>> asReducedGenericModule = GenericPIDModule.fromSyzygies(free, syzygies);
		for (OrderIdealClass ic : elements) {
			ic.asVector = new Vector<>(asGenericModule.asDiagonalVector(ic.asVector).asList()
					.subList(unoptimizedGenerators.size() - free.dimension(), unoptimizedGenerators.size()));
			ic.asVector = asReducedGenericModule.lift(asReducedGenericModule.reduce(ic.asVector));
		}
		asGenericModule = asReducedGenericModule;
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
	public OrderIdealClass getRandomElement() {
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
	public Iterator<OrderIdealClass> iterator() {
		return new Iterator<>() {
			private Iterator<Mod<Vector<IntE>>> it = asIntModule().iterator();

			@Override
			public boolean hasNext() {
				return it.hasNext();
			}

			@Override
			public OrderIdealClass next() {
				return fromVector(asIntModule().lift(it.next()));
			}
		};
	}

	@Override
	public OrderIdealClass neutral() {
		return neutral;
	}

	@Override
	public OrderIdealClass zero() {
		return neutral();
	}

	@Override
	public OrderIdealClass inverse(OrderIdealClass t) {
		if (!inverseMap.containsKey(t.index)) {
			OrderIdealClass cl = getEmbedding(idealGroup.inverse(t.representative));
			inverseMap.put(t.index, cl);
		}
		return inverseMap.get(t.index);
	}

	@Override
	public OrderIdealClass negative(OrderIdealClass s) {
		return inverse(s);
	}

	@Override
	public OrderIdealClass operate(OrderIdealClass t1, OrderIdealClass t2) {
		if (!this.operationMap.containsKey(t1.index)) {
			this.operationMap.put(t1.index, new TreeMap<>());
		}
		if (!this.operationMap.containsKey(t2.index)) {
			this.operationMap.put(t2.index, new TreeMap<>());
		}
		if (!this.operationMap.get(t1.index).containsKey(t2.index)) {
			OrderIdealClass cl = getEmbedding(idealGroup.operate(t1.representative, t2.representative));
			this.operationMap.get(t1.index).put(t2.index, cl);
			this.operationMap.get(t2.index).put(t1.index, cl);
		}
		return operationMap.get(t1.index).get(t2.index);
	}

	@Override
	public OrderIdealClass add(OrderIdealClass s1, OrderIdealClass s2) {
		return operate(s1, s2);
	}

	@Override
	public OrderIdealClass scalarMultiply(IntE t, OrderIdealClass s) {
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
	public boolean isLinearIndependent(List<OrderIdealClass> s) {
		return s.isEmpty();
	}

	private Mod<Vector<IntE>> asIntModuleElement(OrderIdealClass s) {
		return asIntModule().reduce(asVector(s));
	}

	private List<Mod<Vector<IntE>>> asIntModuleElementList(List<OrderIdealClass> s) {
		List<Mod<Vector<IntE>>> result = new ArrayList<>();
		for (OrderIdealClass c : s) {
			result.add(asIntModuleElement(c));
		}
		return result;
	}

	@Override
	public boolean isGeneratingModule(List<OrderIdealClass> s) {
		return asGenericModule.isGeneratingModule(asIntModuleElementList(s));
	}

	@Override
	public List<Vector<IntE>> nonTrivialCombinations(List<OrderIdealClass> s) {
		return asIntModule().nonTrivialCombinations(asIntModuleElementList(s));
	}

	private List<OrderIdealClass> unoptimizedGenerators() {
		Integers z = Integers.z();
		List<OrderIdealClass> generators = new ArrayList<>();
		Reals r = order.getVectorSpace().getValueField().getReals();
		Real sqrtDiscriminantRatio = r.positiveSqrt(
				r.abs(r.divide(r.getInteger(order.discriminant()), r.getInteger(order.numberField().discriminant()))));
		IntE minkowskiLimit = r.multiply(sqrtDiscriminantRatio, order.numberField().minkowskiBound()).roundDown();
		IntE conductor = order.integerConductor();
		NumberFieldIntegers maximalOrder = order.numberField().maximalOrder();
		for (IntE prime : z.setOfPrimes()) {
			if (prime.compareTo(minkowskiLimit) > 0) {
				break;
			}
			if (z.isDivisible(conductor, prime)) {
				continue;
			}
			for (NumberFieldIdeal ideal : maximalOrder.idealsOver(prime)) {
				NumberFieldOrderIdeal restricted = order.restrictFromMaximalOrder(ideal);
				OrderIdealClass generator = getEmbedding(restricted);
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
	public List<OrderIdealClass> getModuleGenerators() {
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
	public List<Vector<IntE>> getSyzygies() {
		return asIntModule().getSyzygies();
	}

	@Override
	public Vector<IntE> asVector(OrderIdealClass s) {
		computeAllElements();
		return s.asVector;
	}
}
