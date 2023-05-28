---
title: "BMB peeves"
---

Here are some of  my opinions about writing.[^1]

## General points

Write plainly and straightforwardly, using concrete words whenever you can. From Robertson Davies' novel *The Manticore* (1972):

> Ramsay always insisted that there was nothing that could not be expressed in the Plain Style if you knew what you were talking about. Everything else was Baroque Style, which he said was not for most people, or Jargon, which was the Devil's work.

I don't have a favorite book on writing. I grew up with Strunk and White's [*Elements of Style*](https://en.wikipedia.org/wiki/The_Elements_of_Style). Although their grammar advice is routinely (and justifiably) ridiculed by the linguists at [Language Log](https://languagelog.ldc.upenn.edu/nll/?p=15509), their stylistic advice is a good starting point &mdash; not *always* to be followed, but worth keeping in mind. I also like Orwell's [Politics and the English Language](https://www.orwellfoundation.com/the-orwell-foundation/orwell/essays-and-other-works/politics-and-the-english-language/) &mdash; also [denigrated by Language Log](https://languagelog.ldc.upenn.edu/nll/?p=992), also not to be taken as gospel, but a useful and enjoyable read.

Prefer the active to the passive voice, and language that attributes agency to vague language. (You might prefer to use the passive voice in a methods section, because you may want to write in a way that downplays *who* followed the procedures; I still prefer the active voice here.) For example, replace "It has been observed that ducks like water (Schmoo et al 2015)" with "Schmoo et al (2015) observed that ducks like water" (but [see below](#ref2))."

Using the first-person singular in scientific writing is controversial, but you should do it if your audience will let you. If you're writing with co-authors you can use "we"; you may also be able to use "we" in the sense of "I and you, the reader" (this approach is common in mathematics textbooks: "Therefore we can easily conclude that ...").

Scientific readers are generally more interested in the current state of the field than in its history. While you may sometimes want to illustrate a controversy or the outline the historical development of some research area, it's generally better to write about what's known *now*, i.e. where your study is starting from. <a name="ref2"></a> Similarly, unless you want to refer to a particular study several times, "Ducks like water (Schmoo et al 2015)" is better than "Schmoo et al (2015) observed that ducks like water".[^2]

When referring to observations or conclusions based on your own work described in the paper, you should say "Ducks like water" rather than "We observed (concluded that, found that, etc.) ducks like water"; it's more succinct, and your readers will be able to infer that such statements are based on your observations and/or logic.[^3].

Don’t introduce your paper by saying that many people have long been interested in the topic: "your paper should introduce the biological topic and explain why it’s interesting and important, not say that other people think the topic is interesting and important" ([Jeremy Fox](https://dynamicecology.wordpress.com/2015/04/20/dont-introduce-your-paper-by-talking-about-how-lots-of-people-are-interested-in-the-topic/)).

Be deliberate about where you put the strong point or punch line of a sentence or paragraph. Depending on the flow of your argument it should go either at the beginning or at the end, rarely in the middle.

Verbs are usually punchier than nouns, including [gerunds](https://en.wikipedia.org/wiki/Gerund#Examples_of_use).

As Strunk and White say, "omit needless words". After writing a draft, go through and see what qualifications and waffles you can delete without changing the meaning of your sentences.

Reading your own prose out loud is a good way to proofread and find awkward phrasing.

## More advice

* Phrase statements positively rather than negatively: "ducks are likely to be larger in the tropics" rather than "ducks are unlikely to be smaller in the tropics", or "our attempt to catch ducks failed" rather than "our attempt to catch ducks did not succeed". (Orwell specifically recommends avoiding the "not un-" construction.)
* I had it beaten into me in high school that "this" should never be used without an antecedent ("Ducks are bigger in the tropics. This means that ...") "This" what? Although I no longer believe it's a hard-and-fast rule, include an antecedent whenever you can ("This enlargement means that ..."); it will make your writing clearer.
* Avoid redundancy, particularly when hedging. Instead of saying "in most cases, ducks tend to be larger" say *either* "in most cases, ducks are larger" *or* "ducks tend to be larger".
* The verb "to be" is vague and can usually be replaced by more concrete and specific phrasing. "There is/are" is a vague way to start a sentence. (Think about *why* you're giving the reader this piece of information, and try to work that motivation into your transition.)
* "May" or "may be" can replace "could suggest that" or "is probably" or "is likely" and other waffle-phrases
* "and" splices: "and" is often used to connect clauses that are only loosely connected: "Ducks are large and when they are approached by humans they make noise". These passages can be stronger if they are split into separate sentences ("Ducks are large. When they are approached by humans they make noise.") or separated with a semicolon ("Ducks are large; when they are approached by humans they make noise.")
* "effect" is overused and vague. "has a positive effect on", "has a positive impact of", "positively influences" → "increases".  "Affect" is better than "effect" (verbs stronger than nouns), but even here you should try to use a more specific word ("increases", "decreases"); "increases" is a lot better than "affects positively". I don't like "impacts" as a verb (I recognize that this preference is old-fashioned/peevish). "Negatively impacts" should be "hurts" or "decreases" or ...
* For ecologists: "density-dependence" is vague. I try to avoid *reifying* density dependence (i.e., convert it from a label that describes a class of phenomena to a *thing*). "Density-dependent mortality" (or fecundity or whatever) is better than "density dependence" alone; the more specific you can be (e.g. attributing the effects to competition for limited resources or nest sites or whatever), the better.

    A similar example from evolutionary biology: "species X has larger female-biased sexual size dimorphism" → "species X has larger females relative to males" or "females are larger in species X"
* Avoid constructions of the type "A and B have effects C and D, respectively" &mdash; you save a few words this way but the reader has to work harder to associate A with C and B with D

## Punctuation

* I like punctuation that isn't part of a quoted phrase or sentence to be *outside* the quotation marks, even though that contradicts most American English style guides.
* I like to denote the introduction/definition of new technical terms by *italicizing them*. "Large ducks, or *gigantaducks*, are found in the tropics."
* Don't overuse italics for emphasis. Italic overuse is one of my own weaknesses; I remove italics from my paper when revising unless they're absolutely necessary.
* Avoid quotation marks except for actual, literal quotations. So-called "scare quotes" are hard to interpret, because they usually mean that you indicate that something shouldn't be taken literally, but it's not precise why or how. Is "duck" a new term that's being defined? Or are you saying "sort of a duck, but not really"? Or ... ?
* Don't use hyphens after adverbs ending with "ly": "partially closed" rather than "partially-closed". (Unless your supervisor/style guide insists).)
* If you use LaTeX, plain double quotes (") may get converted to "fancy" quotes (i.e. “fancy”), sometimes incorrectly (“fancy“): use double-single-backquotes `` for open-double-quotes, and either double-single-quotes '' or plain double-quotes for close-double-quotes.

## Proofreading and citation

Always spell-check. Spell-checking won't catch everything, and you may have to skip over a lot of technical vocabulary that's not in the dictionary, but there's no excuse for not doing it.

Everyone overuses their favorite kinds of punctuation or turns of phrase. Learn what yours are, and use your revision process to thin them out. (For example, I deleted many instances of the word "usually" in this document, keeping only the ones I thought were really necessary.)[^4]

If you are using [BibTeX](http://www.bibtex.org/), check for capitalization! Most journal styles automatically set titles in all lower-case, so you have to protect any words that should be Capitalized or CAPITALIZED with curly brackets {} (some reference managers automatically {protect} {every} {word} {in} {the} {title}, but I prefer that only necessary words are protected).

You should have read, or at least looked at, every reference that you cite. Thanks to the internet it is now plausible that you found a copy of Laplace's "Recherches sur le calcul intégral aux différences infiniment petites, et aux différences finies" (1771) online and read it in the original French[^5], but I will definitely challenge you about it during your committee meeting or thesis defence. If you read something in translation, you should cite the translation/translator; if you feel it's important to cite work that's cited in another work (because it is historically important or foundational), state that explicitly (e.g. "Laplace *Recherches sur le calcul ...* [1771], cited in Schmoo et al [2015]").
I prefer author-date style citation. It is less compact than using numbered references, but it immediately lets readers who are familiar with the literature in your field know what study you're citing, without having to skip to the reference section.

Use a bibliographic citation manager. I prefer [Zotero](https://www.zotero.org/) because it is free, powerful, and convenient, but you can use anything that works for you and your co-authors. Using a citation manager means you can keep track of everything you've read whether or not you use it in your paper; any modern citation manager + document preparation system will make the process of switching reference formats nearly painless.

When citing scholarly work that is available online (e.g. a peer-reviewed journal article or conference proceedings, or an article in a preprint archive), give the full citation information, not just the URL or [DOI](https://www.doi.org/); it's more formal, and readers who are familiar with the literature in your area will be able to know what you're citing without having to click through a link. Do include the URL or DOI in the citation, for convenience. Use the DOI in place of the URL of the publisher's web site wherever possible. 

## Vocabulary/words to avoid

- "utili[sz]e" (use "use" instead)
- "hypothesize" (Can be useful if you are laying out a formal hypothesis, but why not just "think" or "suggest"?). "It has been hypothesized that" is both passive and (IMO) unnecessarily formal. ("Schmoo et al. (2024) hypothesized that" is OK, but "proposed" or "suggested" could be good, less jargony alternatives)
- vague qualifiers: "quite", "very", "extremely". Could you omit them without changing the meaning of your sentence?
- horn-blowing: "importantly", "interestingly". Ideally, your writing will be spare enough that you don't need to emphasize particular results.
- Ditto for "Note that" and "Notably".
- "incredibly" (this is admittedly a peeve/niche opinion, but I can't help reading this word with some of its etymological connotation of "unbelievably"; in any case, it's more informal than I like for scientific writing)
- "Additionally": this sounds to me like you are making a list and don't know how to motivate the later items more specifically than "and I have another thing to tell you"
- avoid contractions ("don't", "won't") in formal writing. 
- "creatures": this word connotes things that were created, which you might want to avoid in scientific writing. Just use "organisms" (or something general enough for your taxon of choice)
- "unique", "complex", "dynamic", "nonlinear": are you using those words to denote something specific, or are you arm-waving about how cool and interesting your subject is?
 
## Organization

* I often write the introduction last, and almost always write the abstract last.
* I hate writing formal outlines, but it works for  some people. I sometimes write a first draft and *then* write an outline by reducing the draft I wrote to a skeleton of bullet points, rearranging the information I put in the first draft into a better order.
* If you find yourself using vague connectors like "furthermore", "additionally", or "in terms of", or if the topic sentences of your paragraphs are vague, you should think more carefully about *why* you are introducing particular facts in a particular order. Can you use more specific language to bring the reader along to the next piece of information?
* If each section in a literature review covers a single reference source you may need to synthesize more, with each section covering several related references that cover the same topic. Otherwise your paper may be repetitive (because multiple references are likely to cover similar topics) and feel list-y.
* Your introduction should include only the key information necessary to interest the reader in your topic and motivate them to read the rest of the paper. This information can include some brief background and context, but more in-depth discussion of the context and comparison to other work should wait until the discussion.
* Most or all of what you say in the discussion should be based on new information that the reader has gained by reading the rest of the paper. Don't repeat information from the introduction unnecessarily.
* It's important to be transparent about the limitations of your study. I usually include these caveats as the second-to-last paragraph of the discussion (followed, of course, by a final paragraph that is upbeat about the paper's value). Try to find a happy medium between dismissive of your work's limitations and over-apologizing for them.

## Caveat

As George Orwell says, "Break any of these rules sooner than say anything outright barbarous."

## [Muphry's law](https://en.wikipedia.org/wiki/Muphry%27s_law)

I'm sure that I have violated some of my own rules while writing this document.

---

Last updated 2023-05-27

[^1]: I've titled the document "peeves" in recognition of the fact that opinions about language and writing are necessarily subjective; the linguists at Language Log, who are also referred to elsewhere in this document, call [prescriptivist](https://dictionary.cambridge.org/dictionary/english/prescriptivist) opinions about language ["peeves"](https://languagelog.ldc.upenn.edu/nll/?p=3144)

[^2]: My mother, Joan Bolker, wrote a book on [*Writing Your Dissertation in Fifteen Minutes a Day*](https://www.amazon.ca/Writing-Your-Dissertation-Fifteen-Minutes/dp/080504891X), edited an anthology of essays about writing called [*The Writer's Home Companion*](https://www.amazon.ca/Writers-Home-Companion-Anthology-Writing/dp/0805048936), and gave advice for medical researchers and writers in [*Writing Medicine*](https://www.lulu.com/en/ca/shop/joan-bolker/writing-medicine/paperback/product-1jz9mq6d.html?page=1&pageSize=4). All are highly recommended, but they focus on the *process* of writing rather than on mechanics like word choice and structure.)

[^3]: This advice contradicts the previous advice about agency, attributing ideas to particular actors. "A foolish consistency is the hobgoblin of little minds" (Emerson: [Wikipedia points out that this is a misuse of the original quotation](https://en.wikipedia.org/wiki/Wikipedia:Emerson_and_Wilde_on_consistency), but I like it anyway.) Or: "Do I contradict myself? Very well then I contradict myself, (I am large, I contain multitudes.)" [Whitman, *Song of Myself*](https://poets.org/poem/song-myself-51)

[^4]: I also like footnotes, but they're rarely used, and rarely a good idea, in scientific writing. (They're common in the humanties.)

[^5]: published in 1776,




