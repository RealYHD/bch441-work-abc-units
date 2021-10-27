# ==   1.3  Task: submit for credit (part 1/2)  ================================
# == Submission - Code to add another philosopher to the datamodel:

pID <- autoincrement(philDB$person)
immanuelKant <- data.frame(id = pID,
                           name = "Immanuel Kant",
                           born = "1724",
                           died = "1804",
                           school = "Enlightenment Philosophy")
philDB$person <- rbind(philDB$person, immanuelKant)

bID = autoincrement(philDB$books)
immanuelKantWork <- data.frame(id = bID,
                               title = "Critique of Pure Reason",
                               published = "1781")
philDB$books <- rbind(philDB$books, immanuelKantWork)
philDB$works <- rbind(philDB$works, data.frame(id = autoincrement(philDB$works), personID = pID, bookID = bID))

bID = autoincrement(philDB$books)
immanuelKantWork <- data.frame(id = bID,
                               title = "Critique of Judgement",
                               published = "1790")
philDB$books <- rbind(philDB$books, immanuelKantWork)
philDB$works <- rbind(philDB$works, data.frame(id = autoincrement(philDB$works), personID = pID, bookID = bID))

# == Submission: Code to list the philosophical schools in alphabetical order as well as their respective books in alphabetical order.

schools <- unique(philDB$person$school)
schools <- sort(schools)

for (s in schools) {
  cat(sprintf("%s\n", s))
  authors = which(philDB$person$school == s)
  for (author in authors) {
    works = which(philDB$works$personID == author)
    for (work in works) {
      bookId = which(philDB$books$id == philDB$works$bookID[work])
      cat(sprintf("\t%s - (%s)\n", philDB$books$title[bookId], philDB$books$published[bookId]))
    }
  }
}