package util;

import java.lang.reflect.Array;
import java.util.Collection;
import java.util.Comparator;
import java.util.Iterator;
import java.util.SortedMap;
import java.util.SortedSet;

public class SortedKeySet<K, V> implements SortedSet<K> {
	private SortedMap<K, V> map;
	
	public SortedKeySet(SortedMap<K, V> map) {
		this.map = map;
	}

	@Override
	public int size() {
		return map.size();
	}

	@Override
	public boolean isEmpty() {
		return map.isEmpty();
	}

	@Override
	public boolean contains(Object o) {
		return map.containsKey(o);
	}

	@Override
	public Iterator<K> iterator() {
		return map.keySet().iterator();
	}

	@Override
	public Object[] toArray() {
		Object[] array = new Object[map.size()];
		Iterator<K> it = iterator();
		for (int i = 0; i < map.size(); i++) {
			array[i] = it.next();
		}
		return array;
	}

	@SuppressWarnings("unchecked")
	@Override
	public <T> T[] toArray(T[] a) {
		if (a.length <= map.size()) {
			a = (T[])Array.newInstance(a.getClass().arrayType(), map.size());
		}
		Iterator<K> it = iterator();
		for (int i = 0; i < map.size(); i++) {
			a[i] = (T)it.next();
		}
		if (a.length > map.size()) {
			a[map.size()] = null;
		}
		return a;
	}

	@Override
	public boolean add(K e) {
		throw new UnsupportedOperationException("Unmodifcable Set!");
	}

	@Override
	public boolean remove(Object o) {
		throw new UnsupportedOperationException("Unmodifcable Set!");
	}

	@Override
	public boolean containsAll(Collection<?> c) {
		for (Object t : c) {
			if (!contains(t)) {
				return false;
			}
		}
		return true;
	}

	@Override
	public boolean addAll(Collection<? extends K> c) {
		throw new UnsupportedOperationException("Unmodifcable Set!");
	}

	@Override
	public boolean retainAll(Collection<?> c) {
		throw new UnsupportedOperationException("Unmodifcable Set!");
	}

	@Override
	public boolean removeAll(Collection<?> c) {
		throw new UnsupportedOperationException("Unmodifcable Set!");
	}

	@Override
	public void clear() {
		throw new UnsupportedOperationException("Unmodifcable Set!");
	}

	@Override
	public Comparator<? super K> comparator() {
		return map.comparator();
	}

	@Override
	public SortedSet<K> subSet(K fromElement, K toElement) {
		return new SortedKeySet<K, V>(map.subMap(fromElement, toElement));
	}

	@Override
	public SortedSet<K> headSet(K toElement) {
		return new SortedKeySet<K, V>(map.headMap(toElement));
	}

	@Override
	public SortedSet<K> tailSet(K fromElement) {
		return new SortedKeySet<K, V>(map.tailMap(fromElement));
	}

	@Override
	public K first() {
		return map.firstKey();
	}

	@Override
	public K last() {
		return map.lastKey();
	}
}
