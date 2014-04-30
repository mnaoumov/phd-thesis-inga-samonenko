using System.Web.Mvc;

namespace PhdThesisIngaSamonenko.Controllers
{
    public class KMDcodeController : Controller
    {
        public ActionResult Index()
        {
            return File("KMD.R", "text/plain");
        }

    }
}
