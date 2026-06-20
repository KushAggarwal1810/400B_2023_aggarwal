const startScreen = document.querySelector("#startScreen");
const startButton = document.querySelector("#startButton");
const song = document.querySelector("#song");
const musicButton = document.querySelector("#musicButton");
const musicText = document.querySelector("#musicText");
const featuredPhoto = document.querySelector("#featuredPhoto");
const photoCount = document.querySelector("#photoCount");
const photoTitle = document.querySelector("#photoTitle");
const thumbGrid = document.querySelector("#thumbGrid");
const prevPhoto = document.querySelector("#prevPhoto");
const nextPhoto = document.querySelector("#nextPhoto");

const videos = document.querySelectorAll("video");

const photos = Array.from({ length: 28 }, (_, index) => {
  const number = String(index + 1).padStart(2, "0");
  return `assets/site/photo-${number}.jpeg`;
});

const captions = [
  "A memory that stays close.",
  "A smile that makes home feel complete.",
  "A moment held safely forever.",
  "The kind of love that never asks for credit.",
  "Papa, your presence changes the whole room.",
  "A picture, a heartbeat, a blessing.",
  "Every ordinary day became special with you.",
  "Family feels stronger because of you.",
  "Your love is the quiet center of everything.",
  "The best memories always have you in them.",
  "A forever kind of photograph.",
  "A chapter I never want to forget.",
  "Thank you for every unseen sacrifice.",
  "The warmth we keep coming back to.",
  "Proof that love can be gentle and strong.",
  "A moment full of pride and prayer.",
  "The face of home.",
  "Your blessings travel with us everywhere.",
  "A still frame of a lifelong bond.",
  "Love, respect, and gratitude in one photo.",
  "A memory made brighter by you.",
  "The person behind so many smiles.",
  "Papa, you are our safest place.",
  "A photo that says what words cannot.",
  "The heart of our family.",
  "A blessing we get to call Papa.",
  "Forever grateful, forever proud.",
  "Papa, even the flowers wanted you in the frame."
];

let currentPhoto = 0;

function setMusicState(isPlaying) {
  document.body.classList.toggle("music-paused", !isPlaying);
  musicText.textContent = isPlaying ? "Music on" : "Music off";
}

async function playSong() {
  try {
    await song.play();
    setMusicState(true);
  } catch {
    setMusicState(false);
  }
}

function showPhoto(index) {
  currentPhoto = (index + photos.length) % photos.length;
  featuredPhoto.classList.add("is-changing");

  window.setTimeout(() => {
    featuredPhoto.src = photos[currentPhoto];
    photoCount.textContent = `${currentPhoto + 1} / ${photos.length}`;
    photoTitle.textContent = captions[currentPhoto];
    featuredPhoto.classList.remove("is-changing");
  }, 160);

  thumbGrid.querySelectorAll("button").forEach((button, buttonIndex) => {
    button.classList.toggle("active", buttonIndex === currentPhoto);
  });
}

function buildThumbnails() {
  photos.forEach((src, index) => {
    const button = document.createElement("button");
    button.type = "button";
    button.setAttribute("aria-label", `Show memory ${index + 1}`);

    const img = document.createElement("img");
    img.src = src;
    img.alt = "";

    button.appendChild(img);
    button.addEventListener("click", () => showPhoto(index));
    thumbGrid.appendChild(button);
  });
}

startButton.addEventListener("click", async () => {
  window.scrollTo({ top: 0, left: 0, behavior: "auto" });
  startScreen.classList.add("hidden");
  await playSong();
});

musicButton.addEventListener("click", async () => {
  if (song.paused) {
    await playSong();
  } else {
    song.pause();
    setMusicState(false);
  }
});

videos.forEach((video) => {
  function resumeSongAfterVideo() {
    const anotherVideoPlaying = Array.from(videos).some((item) => item !== video && !item.paused && !item.ended);

    if (!anotherVideoPlaying && song.dataset.pausedForVideo === "true") {
      song.dataset.pausedForVideo = "";
      playSong();
    }
  }

  video.addEventListener("play", () => {
    if (!song.paused) {
      song.dataset.pausedForVideo = "true";
      song.pause();
      setMusicState(false);
    }

    videos.forEach((otherVideo) => {
      if (otherVideo !== video) {
        otherVideo.pause();
      }
    });
  });

  video.addEventListener("pause", resumeSongAfterVideo);
  video.addEventListener("ended", resumeSongAfterVideo);
});

prevPhoto.addEventListener("click", () => showPhoto(currentPhoto - 1));
nextPhoto.addEventListener("click", () => showPhoto(currentPhoto + 1));

document.addEventListener("visibilitychange", () => {
  if (document.hidden && !song.paused) {
    song.dataset.wasPlaying = "true";
    song.pause();
    setMusicState(false);
  } else if (!document.hidden && song.dataset.wasPlaying === "true") {
    song.dataset.wasPlaying = "";
    playSong();
  }
});

buildThumbnails();
showPhoto(0);
setMusicState(false);
