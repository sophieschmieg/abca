package util;

import java.util.Optional;
import java.util.SortedMap;
import java.util.TreeMap;

public class Cache<K extends Comparable<? super K>, V> {
	private Entry first;
	private Entry last;
	private int size;
	private int capacity;
	private SortedMap<K, Entry> index;

	public Cache(int capacity) {
		this.capacity = capacity;
		this.index = new TreeMap<>();
		this.size = 0;
		this.first = null;
		this.last = null;
	}

	public void insert(K key, V value) {
		remove(key);
		size++;
		while (size > capacity) {
			remove(last.key);
		}
		Entry newEntry = new Entry();
		newEntry.key = key;
		newEntry.value = value;
		index.put(key, newEntry);
		if (first == null) {
			first = newEntry;
			last = newEntry;
			return;
		}
		first.prev = newEntry;
		newEntry.next = first;
		first = newEntry;
	}

	public Optional<V> lookup(K key) {
		Entry entry = index.get(key);
		if (entry == null) {
			return Optional.empty();
		}
		if (entry != first) {
			if (entry == last) {
				last = entry.prev;
			} else {
				entry.next.prev = entry.prev;
			}
			first.prev = entry;
			entry.prev.next = entry.next;
			entry.next = first;
			entry.prev = null;
			first = entry;
		}
		return Optional.of(entry.value);
	}

	public void remove(K key) {
		Entry entry = index.remove(key);
		if (entry == null) {
			return;
		}
		if (entry == first) {
			first = entry.next;
		} else {
			entry.prev.next = entry.next;
		}
		if (entry == last) {
			last = entry.prev;
		} else {
			entry.next.prev = entry.prev;
		}
		size--;
	}

	public int size() {
		return size;
	}

	private class Entry {
		private K key;
		private V value;
		private Entry prev;
		private Entry next;
	}
}
