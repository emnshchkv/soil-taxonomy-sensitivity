# Поиск по своей библиотеке: taxonomic diversity of soil bacterial communities

[источник] OpenAlex, запрос: https://api.openalex.org/works?search=taxonomic+diversity+of+soil+bacterial+communities&per_page=50&filter=has_abstract%3Atrue&select=id%2Cdoi%2Cdisplay_name%2Cpublication_year%2Ccited_by_count%2Cabstract_inverted_index, дата: 2026-10-18, работ: 20

## dependence of bacterial composition on soil type

| место | BM25 | вектор |
|---|---|---|
| 1 | Contrasting Soil Bacterial Community, Diversity, and Functio (W2886861276) | Soil textural heterogeneity impacts bacterial but not fungal (W3009530196) |
| 2 | Soil textural heterogeneity impacts bacterial but not fungal (W3009530196) | Total bacterial diversity in soil and sediment communities—A (W2018218203) |
| 3 | Total bacterial diversity in soil and sediment communities—A (W2018218203) | Pyrosequencing-Based Assessment of Soil pH as a Predictor of (W2127602818) |

[написано] Что подошло и чем выдачи отличаются: по теме подходят работы про текстуру почвы (W3009530196) и рН как предиктор состава сообщества (W2127602818). Это прямые примеры зависимости бактериального состава от свойств почвы. Выдачи совпадают на двух из трех позиций. Вектор поднял работу про текстуру на первое место и вместо общей статьи «Contrasting Soil Bacterial Community…» дал работу про рН, которая точнее отвечает на вопрос о факторах среды. BM25 ставит на первое место статью с самым частым совпадением слов, а по заголовку не видно, что она именно про тип почвы.

## species diversity of bacterial communities in Moscow

| место | BM25 | вектор |
|---|---|---|
| 1 | The Diversity of Bacterial Communities in Urban Soils (W2896176982) | The Diversity of Bacterial Communities in Urban Soils (W2896176982) |
| 2 | Effects of plant community history, soil legacy and plant di (W3165382877) | Diversity and Taxonomy of Soil Bacterial Communities in Urba (W4308440034) |
| 3 | Using ecological diversity measures with bacterial communiti (W2129297365) | Total bacterial diversity in soil and sediment communities—A (W2018218203) |

[написано] Что подошло и чем выдачи отличаются: обе выдачи ставят на первое место одну и ту же работу про городские почвы (W2896176982), и это самое близкое к запросу. Заголовок не показывает, идет ли речь именно о Москве, это надо проверить по аннотации. Вторая позиция у вектора, W4308440034 (разнообразие и таксономия бактерий городских почв), тоже по теме. У BM25 вторая и третья позиции совпали по словам «community», «diversity» и «bacterial». Работа про растительные сообщества (W3165382877) и методическая работа про меры разнообразия (W2129297365) к городской почве отношения почти не имеют. Здесь вектор заметно лучше.


## classifiers of bacterial communities

| место | BM25 | вектор |
|---|---|---|
| 1 | Functional traits dominate the diversity-related selection o (W2501364557) | Taxonomic and functional diversity of atrazine‐degrading bac (W2141155745) |
| 2 | Taxonomic and functional diversity of atrazine‐degrading bac (W2141155745) | Pyrosequencing-Based Assessment of Soil pH as a Predictor of (W2127602818) |
| 3 | Pyrosequencing-Based Assessment of Soil pH as a Predictor of (W2127602818) | Contrasting Soil Bacterial Community, Diversity, and Functio (W2886861276) |

[написано] Что подошло и чем выдачи отличаются: по заголовкам ни одна работа не посвящена классификаторам в смысле алгоритмов таксономического отнесения. Все результаты про разнообразие и функциональные признаки сообществ, так что попадание по сути отсутствует. Запрос двусмысленный: «classifiers» может означать и ML-классификаторы, и признаки, по которым сообщества делят на типы. Выдачи пересекаются на двух работах (W2141155745, W2127602818), и обе системы вытягивают общую тематическую группу, а не ответ на вопрос. Вероятнее всего, в библиотеке такой работы просто нет.


## differences between naive Bayes classifiers for bacteria

| место | BM25 | вектор |
|---|---|---|
| 1 | Contrasting Soil Bacterial Community, Diversity, and Functio (W2886861276) | Diversity and Taxonomy of Soil Bacterial Communities in Urba (W4308440034) |
| 2 | Taxonomic Diversity and Functional Traits of Soil Bacterial  (W4393899294) | Contrasting Soil Bacterial Community, Diversity, and Functio (W2886861276) |
| 3 | Functional traits dominate the diversity-related selection o (W2501364557) | Pyrosequencing-Based Assessment of Soil pH as a Predictor of (W2127602818) |

[написано] Что подошло и чем выдачи отличаются: ничего не подошло, потому что работ про наивный байесовский классификатор в библиотеке нет. BM25 зацепился за слова «bacterial» и «taxonomic». Вектор выше ставит работу про таксономию (W4308440034), то есть ближе к теме отнесения к таксонам, но это сходство тематическое, а не по методу. Для содержательного ответа нужна другая выборка, например по запросу «naive Bayes classifier 16S rRNA taxonomy assignment» (RDP, q2-feature-classifier).


## symbiosis and competition among soil bacteria

| место | BM25 | вектор |
|---|---|---|
| 1 | Contrasting Soil Bacterial Community, Diversity, and Functio (W2886861276) | Contrasting Soil Bacterial Community, Diversity, and Functio (W2886861276) |
| 2 | Taxonomic Diversity and Functional Traits of Soil Bacterial  (W4393899294) | The Diversity of Bacterial Communities in Urban Soils (W2896176982) |
| 3 | Functional traits dominate the diversity-related selection o (W2501364557) | Effects of plant community history, soil legacy and plant di (W3165382877) |

[написано] Что подошло и чем выдачи отличаются: прямых работ про симбиоз и конкуренцию по заголовкам нет. Обе выдачи совпадают на первой позиции (W2886861276), но это общая работа о сообществе и функциях, и взаимодействия там могут быть только побочной темой. BM25 дальше выдает работы про функциональные признаки, которые теоретически могут касаться отбора и конкуренции (W4393899294, W2501364557). Вектор подмешивает работу про городские почвы, она по теме не подходит, и работу про наследие растительного сообщества (W3165382877), которая ближе всего к межвидовым взаимодействиям, пусть и через растения. Что именно подходит, можно решить только по аннотациям.
