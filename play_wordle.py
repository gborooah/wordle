import sys
from collections import Counter
from statistics import mean


def solve(word):
    """
    Solve one game, outputting guesses and clues.
    :param word: solution
    :return:
    """
    valid_solutions = set(load_solutions())
    if word not in valid_solutions:
        print(word)
        raise ValueError("Invalid solution")
    game = Game(word)
    solver = Solver(max_turns=game.max_turns, first_guess="raise", second_guesses=load_second_guesses())
    game.autoplay(solver)


def solve_all(first_guess="raise"):
    """
    Solve all games, outputting solution and number of turns needed.
    Reads input data once for efficiency.
    :return:
    """
    solutions = load_solutions()
    guesses = load_guesses()
    second_guesses = load_second_guesses()
    for s in solutions:
        print(s, file=sys.stderr)
        game = Game(s, debug=0)
        solver = Solver(max_turns=game.max_turns, first_guess=first_guess, solutions=solutions, guesses=guesses,
                        second_guesses=second_guesses)
        game.autoplay(solver)


def generate_second_guess(first_guess_str):
    """
    Generate best second guesses for each possible clue from given first guess.
    :param first_guess_str:
    :return:
    """
    solutions = load_solutions()
    first_guess = Word(first_guess_str)
    clues_seen = set()
    for s in solutions:
        print(s, file=sys.stderr)
        solver = Solver(max_turns=3)
        first_clue = Clue(first_guess, s)
        if first_clue in clues_seen:
            continue
        clues_seen.add(first_clue)
        solver.add_clue(first_clue)
        second_guess = solver.generate_guess(2)
        print(first_clue.to_data_str() + " " + str(second_guess))


def find_worst_clues():
    guesses = load_guesses()
    solver = Solver()
    for g in guesses:
        solver.assess_guess(g, debug=2)


class Word:
    def __init__(self, word_str):
        self.word = list(word_str.upper())
        if len(self.word) != 5:
            raise ValueError(word_str)
        self.word_bag = Counter()
        for c in self.word:
            self.word_bag[c] += 1

    def __str__(self):
        return "".join([c.upper() for c in self.word])

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        return str(self) == str(other)


class Clue:
    # returns parsed clue -1 unset, 0 miss, 1 partial match, 2 exact match
    def __init__(self, guess: Word, solution: Word = None, pattern=None):
        self.guess = guess
        if solution is None:
            self.data = zip(guess.word, map(lambda x: int(x), list(pattern)))
            return
        self.data = [(c, -1) for c in guess.word]
        seen = Counter()
        for (idx, pair) in enumerate(zip(guess.word, solution.word)):
            c = pair[0]
            d = pair[1]
            if c == d:
                self.data[idx] = (c, 2)
                seen[c] += 1
            else:
                self.data[idx] = (c, 0)
        for (idx, pair) in enumerate(zip(guess.word, solution.word)):
            c = pair[0]
            d = pair[1]
            if c != d:
                if solution.word_bag[c] - seen[c] > 0:
                    self.data[idx] = (c, 1)
                    seen[c] += 1

    # tests if word satisfies clue
    # not what's needed, need to know if this possible would give this clue
    def match(self, word_obj: Word):
        seen = Counter()
        for idx in range(5):
            if self.data[idx][1] == 2:
                if self.data[idx][0] != word_obj.word[idx]:
                    return False
                else:
                    seen[word_obj.word[idx]] += 1
        for idx in range(5):
            if self.data[idx][1] == 1:
                ch = self.data[idx][0]
                if word_obj.word[idx] == ch or word_obj.word_bag[ch] - seen[ch] <= 0:
                    return False
                else:
                    seen[ch] += 1
        for idx in range(5):
            if self.data[idx][1] == 0:
                ch = self.data[idx][0]
                if word_obj.word[idx] == ch or word_obj.word_bag[ch] > seen[ch]:
                    return False
        return True

    def matches(self, possible: Word):
        possible_clue = Clue(self.guess, possible)
        return self == possible_clue

    def exact_match(self):
        return all([t[1] == 2 for t in self.data])

    def __str__(self):
        def pr(tup):
            if tup[1] == 0:
                return '\033[2;37;41m' + tup[0].upper() + '\033[0;0m'
            elif tup[1] == 1:
                return '\033[2;37;44m' + tup[0].upper() + '\033[0;0m'
            else:
                return '\033[2;37;42m' + tup[0].upper() + '\033[0;0m'

        output_list = [pr(d) for d in self.data]
        return "".join(output_list)

    def to_data_str(self):
        return str(self.guess) + " " + "".join([str(t[1]) for t in self.data])

    def __hash__(self):
        return hash("".join(t[0] + str(t[1]) for t in self.data))

    def __eq__(self, other):
        for s, t in zip(self.data, other.data):
            if s[0] != t[0] or s[1] != t[1]:
                return False
        return True


class Solver:

    def __init__(self, max_turns=6, first_guess="raise", solutions=None, guesses=None, second_guesses=None):
        """
        Solver for a game
        :param max_turns: maximum number of turns allowed. Important only to identify last turn, on which solver
        makes a guess of possibles, rather than optimize for some future guess
        :param first_guess: optional pre-computed first guess (best guess should be the same for all games)
        :param solutions: list of valid solutions, useful to pass in if solving multiple games.
        :param guesses:  list of valid guesses
        :param second_guesses: pre-computed second guesses, as dict: clue -> guess
        """
        self.last_turn = max_turns
        self.first_guess = None if first_guess is None else Word(first_guess)
        self.possibles = load_solutions() if solutions is None else solutions
        self.common = set(self.possibles)
        self.guesses = load_guesses() if guesses is None else guesses
        self.second_guesses = dict() if second_guesses is None else second_guesses
        self.last_clue = None

    def assess_guess(self, guess: Word, debug=0):
        """
        Find how many possibles the worst-case clue leaves.
        :param guess: guess word to assess
        :param debug: level of output to produce
        :return: worst case number of possibles after this guess; tiebreaker expected number of possibles
        """
        clue_counts = Counter()
        for p in self.possibles:
            clue = Clue(guess, p)
            clue_counts[clue] += 1
        if debug > 1:
            worst_clue = max(clue_counts, key=clue_counts.get)
            print(worst_clue.to_data_str() + " " + str(clue_counts.get(worst_clue)))
        return max(clue_counts.values()), mean(clue_counts.values())

    def generate_guess(self, turn, debug=0):
        """
        Generates best guess,
        that is the guess which gives worst-case best chance of success if next clue is the last.
        :param turn: First two turns can be pre-computed, and last turn is a guess of the possibles
        :param debug: Gives more or less information about solving process.
        :return: word to guess
        """
        if turn == 0 and self.first_guess is not None:
            return self.first_guess
        if turn == 1 and self.last_clue in self.second_guesses:
            return self.second_guesses.get(self.last_clue)
        if turn == self.last_turn:
            # blind guess
            return self.possibles[0]
        if len(self.possibles) == 1:
            return self.possibles[0]
        # main algo. Assess each guess, pick best
        best_score = (10000, 10000)
        best_guess = Word("humph")
        for g in self.guesses:
            score = self.assess_guess(g)
            if score < best_score or (score == best_score and not self.is_common(best_guess) and self.is_common(g)):
                best_guess = g
                best_score = score
                if debug > 1:
                    print(str(best_guess) + ": " + str(score))
        return best_guess

    def add_clue(self, clue: Clue):
        """
        Update possibles given a clue.
        :param clue: clue
        :return:
        """
        self.last_clue = clue
        self.possibles = [p for p in self.possibles if clue.matches(p)]

    def is_common(self, word):
        """
        Is a potential guess also a valid solution. Helps make the solver seem more human.
        :param word: potentially obscure word
        :return: True if in solution set
        """
        return word in self.common


class Game:
    def __init__(self, solution: Word, max_turns=6, debug=1):
        """
        A Wordle game, along with functionality to autoplay against a solver.
        :param solution:
        :param max_turns:
        :param debug:
        """
        self.solution = solution
        self.max_turns = max_turns
        self.debug = debug
        self.turn = 0

    def assess(self, guess: Word):
        return Clue(guess, self.solution)

    def autoplay(self, solver: Solver):
        while self.turn < self.max_turns:
            guess = solver.generate_guess(self.turn)
            clue = self.assess(guess)
            if self.debug > 0:
                print(str(clue) + " " + str(len(solver.possibles)))
            self.turn += 1
            if clue.exact_match():
                # print("you win")
                print(str(self.solution) + " " + str(self.turn))
                return
            solver.add_clue(clue)
        print(str(self.solution) + " FAIL")


def load_words(filename):
    with open(filename) as file:
        return [Word(w.rstrip()) for w in sorted(file.readlines())]


def load_guesses():
    return load_words("resources/guesses.txt")


def load_solutions():
    return load_words("resources/solutions.txt")


def load_second_guesses():
    second_guesses = dict()
    with open("resources/second_guesses.txt") as file:
        for line in file:
            clue_word, pattern, next_guess = line.split()
            clue = Clue(Word(clue_word), pattern=pattern)
            second_guesses[clue] = Word(next_guess)
    return second_guesses


def main():
    # solve(Word(sys.argv[1]))
    # generate_second_guess("raise")
    solve_all()


if __name__ == "__main__":
    main()
