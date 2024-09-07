from spoofbot import Browser, Firefox, Windows
from spoofbot.adapter.file import FileCache

import verdantguardian
from verdantguardian.utils.path import CWD

class VGSession(Browser):
    def user_agent(self) -> str:
        return f"VerdantGuardian {verdantguardian.__version__}"

SESS = VGSession(os=Windows(), adapter=FileCache(CWD / '.cache'))
