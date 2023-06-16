package fields.helper;

import java.math.BigInteger;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

import fields.exceptions.InfinityException;
import fields.helper.PermutationGroup.Permutation;
import fields.interfaces.Group;
import util.MiscAlgorithms;

public class PermutationGroup implements Group<Permutation> {
	public static class Permutation extends AbstractElement<Permutation> {
		private int[] array;

		private Permutation(int[] permutation) {
			this.array = permutation;
		}

		public int[] asArray() {
			return array;
		}

		@Override
		public int compareTo(Permutation o) {
			return Arrays.compare(array, o.array);
		}

		public int evaluate(int t) {
			return array[t];
		}

		public int[] act(int[] input) {
			if (input.length != array.length) {
				throw new ArithmeticException("Unequal length");
			}
			int[] result = new int[input.length];
			for (int i = 0; i < input.length; i++) {
				result[i] = input[evaluate(i)];
			}
			return result;
		}

		public <T> List<T> act(List<T> input) {
			if (input.size() != array.length) {
				throw new ArithmeticException("Unequal length");
			}
			List<T> result = new ArrayList<>();
			for (int i = 0; i < input.size(); i++) {
				result.add(input.get(evaluate(i)));
			}
			return result;
		}
	}

	private Permutation neutral;
	private int size;

	public PermutationGroup(int size) {
		this.size = size;
	}

	@Override
	public Exactness exactness() {
		return Exactness.EXACT;
	}

	public Permutation getPermutation(int[] array) {
		return new Permutation(array);
	}

	public Permutation getPermutation(List<Integer> list) {
		int[] array = new int[size];
		for (int i = 0; i < size; i++) {
			array[i] = list.get(i);
		}
		return new Permutation(array);
	}

	@Override
	public Permutation getRandomElement() {
		Random rng = new SecureRandom();
		int[] result = new int[size];
		for (int i = 0; i < size; i++) {
			result[i] = i;
		}
		for (int i = size - 1; i > 0; i--) {
			int index = rng.nextInt(i + 1);
			int tmp = result[index];
			result[index] = result[i];
			result[i] = tmp;
		}
		return getPermutation(result);
	}

	@Override
	public boolean isFinite() {
		return true;
	}

	private static BigInteger factorial(int n) {
		if (n == 0) {
			return BigInteger.ONE;
		}
		return BigInteger.valueOf(n).multiply(factorial(n - 1));
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		return factorial(size);
	}

	@Override
	public Iterator<Permutation> iterator() {
		List<Integer> list = new ArrayList<>();
		for (int i = 0; i < size; i++) {
			list.add(i);
		}
		return new Iterator<>() {
			private Iterator<List<Integer>> it = MiscAlgorithms.permutations(list).iterator();

			@Override
			public boolean hasNext() {
				return it.hasNext();
			}

			@Override
			public Permutation next() {
				return getPermutation(it.next());
			}
		};
	}

	@Override
	public Permutation neutral() {
		if (neutral == null) {
			int[] result = new int[size];
			for (int i = 0; i < size; i++) {
				result[i] = i;
			}
			neutral = getPermutation(result);
		}
		return neutral;
	}

	@Override
	public Permutation inverse(Permutation t) {
		int[] result = new int[size];
		for (int i = 0; i < size; i++) {
			result[t.array[i]] = i;
		}
		return getPermutation(result);
	}

	@Override
	public Permutation operate(Permutation t1, Permutation t2) {
		int[] result = new int[size];
		for (int i = 0; i < size; i++) {
			result[i] = t1.evaluate(t2.evaluate(i));
		}
		return getPermutation(result);
	}
}
