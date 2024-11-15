$(document).on("submit", "#checkout-form-conj", function (e) {
    return (
        e.preventDefault,
        $.ajax({
            type: "post",
            url: '{% url "conjugation" %}',
            data: $(this).serialize(),
            success: function (e) {
                var t = [e.tense1, e.tense2, e.tense3, e.tense4, e.tense5, e.tense6, e.tense7, e.tense8, e.tense9, e.tense10, e.tense11, e.tense12, e.tense13, e.tense14, e.tense15, e.tense16, e.tense17, e.tense18, e.tense19, e.tensepc],
                    n = [
                        e.present,
                        e.present_cont,
                        e.present_perfect,
                        e.present_perfect_cont,
                        e.past,
                        e.past_perfect,
                        e.past_perfect_cont,
                        e.future,
                        e.future_cont,
                        e.future_perfect,
                        e.future_perfect_cont,
                        e.conditional_present,
                        e.conditional_present_cont,
                        e.conditional_perfect,
                        e.conditional_perfect_cont,
                        e.imperative,
                        e.infinitive,
                        e.participle_present,
                        e.participle_past,
                        e.past_continuous
                    ],
                    i = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T"];
                document.getElementById("infini").innerHTML = e.infinitive[0][1];
                for (var s = 0; s < t.length; s++) document.getElementById("titleconjug" + s).innerHTML = t[s];
                for (var r = 0; r < t.length; r++)
                    for (var c = 0; c < n[r].length; c++)
                        ["P", "Q", "R", "S"].includes(i[r]) || ((document.getElementById(i[r] + "cjtd" + c).innerHTML = n[1][c][0]), $("#" + i[r] + "cjtd" + c).css({ color: "#006ad3" })),
                            3 == n[r][c][1].split(" ").length
                                ? ((document.getElementById(i[r] + "cjtda" + c).innerHTML = n[r][c][1].split(" ")[0]), (document.getElementById(i[r] + "cjtdb" + c).innerHTML = n[r][c][1].split(" ")[1]))
                                : 2 == n[r][c][1].split(" ").length
                                ? (document.getElementById(i[r] + "cjtda" + c).innerHTML = n[r][c][1].split(" ")[0])
                                : 4 == n[r][c][1].split(" ").length
                                ? ((document.getElementById(i[r] + "cjtda" + c).innerHTML = n[r][c][1].split(" ")[0]),
                                  (document.getElementById(i[r] + "cjtdb" + c).innerHTML = n[r][c][1].split(" ")[1]),
                                  (document.getElementById(i[r] + "cjtdc" + c).innerHTML = n[r][c][1].split(" ")[2]))
                                : 5 == n[r][c][1].split(" ").length &&
                                  ((document.getElementById(i[r] + "cjtda" + c).innerHTML = n[r][c][1].split(" ")[0]),
                                  (document.getElementById(i[r] + "cjtdb" + c).innerHTML = n[r][c][1].split(" ")[1]),
                                  (document.getElementById(i[r] + "cjtdc" + c).innerHTML = n[r][c][1].split(" ")[2]),
                                  (document.getElementById(i[r] + "cjtdd" + c).innerHTML = n[r][c][1].split(" ")[3]));
            },
            error: function (e, t, n) {
                console.log(e.status + ": " + e.responseText);
            },
        }),
        !1
    );
});
